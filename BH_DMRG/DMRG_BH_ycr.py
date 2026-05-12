import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt
import multiprocessing as mp
import time



#全局常量
#---------------------------------------------------------------------------------------------------------------
n_max=4  # 单格点占据数截断
d = n_max + 1      #局部物理维数（0,1,2,...,n_max）（每个格点的 Hilbert 空间维数，自旋1/2 对应 2）
D = 625  # 辅助bond 维数上限（初始为 1）
energy_tol=1e-7    #目标基态能量变化
target_trunc=1e-8  #目标svd截断误差
lmd=10       #惩罚参数

#t_list = np.linspace(0.01, 0.25, 15) 
t_list =np.arange(0.02, 0.26, 0.02)
#t_list = np.array([0.02,0.03])
L_list = np.array([32,64,128])
#有限尺寸mu(t,L)
mup_finite_list=np.zeros((len(t_list),len(L_list)))
muh_finite_list=np.zeros((len(t_list),len(L_list)))
#热力学极限下的mu(t)
mup_list=np.zeros_like(t_list)
muh_list=np.zeros_like(t_list)





# MPS（填充n=1，粒子激发，空穴激发）
#---------------------------------------------------------------------------------------------------------------
def generate_MPS(L):    
    # 每个格点张量的形状：(物理维数 d, 左辅助维数, 右辅助维数)，辅助维数随优化过程中s截断动态增加
    # 初始直积态，每个格点都处于 |1⟩（即一个玻色子）

    # 构造单个格点的张量：形状 (d, 1, 1)
    initial_A = np.zeros((d, 1, 1), dtype=float)
    initial_A[1, 0, 0] = 1.0   # 索引 1 对应占据数 |1⟩

    # 构建整个 MPS：所有格点都用同一个张量（直积态，无纠缠）
    MPS = [initial_A] * int(L)
    return MPS

def generate_MPS_p(L): 
    MPS=[]   
    for i in range(L):
        A = np.zeros((d, 1, 1))
        if i == L // 2:        # L // 2 就是 L 除以 2 后向下取整
            A[2, 0, 0] = 1.0   # 系统中央放 |2⟩
        else:
            A[1, 0, 0] = 1.0   # 其余格点 |1⟩
        MPS.append(A)
    return MPS

def generate_MPS_h(L):    
    MPS=[]   
    for i in range(L):
        A = np.zeros((d, 1, 1))
        if i == L // 2:        
            A[0, 0, 0] = 1.0   # 系统中央放 |0⟩
        else:
            A[1, 0, 0] = 1.0   # 其余格点 |1⟩
        MPS.append(A)
    return MPS


# MPO
#---------------------------------------------------------------------------------------------------------------

# 单格点算符
#---------------------------------------------------------------------------------------------------------------
I = np.identity(d)               # 单位算符
Z = np.zeros((d, d))             # 零算符（用于 MPO 占位）

a = np.zeros((d, d))

for i in range(0,n_max): 
    a[i][i+1]=np.sqrt(i+1)
ad =a.transpose()

n= np.zeros((d, d))
for i in range(0,d):
    n[i][i]=i 

h=0.5*(n@(n-I))
sqrtlmdn=np.sqrt(lmd)*n
lmdn2=lmd*(n@n)


# 总粒子数相关的算符
#---------------------------------------------------------------------------------------------------------------
def generate_MPO_N(L):
    W = np.array([[I, n],
                [Z, I]])
    Wfirst = np.array([[I,n]])       # 左边界 MPO（只有第一行）
    #形状 (1, 2, d, d)

    Wlast = np.array([[n],[I]])        # 右边界 MPO（只有最后一列）

    MPO = [Wfirst] + ([W] * (L - 2)) + [Wlast]           # 组装完整 MPO 列表

    return MPO
def generate_MPO_lmdN2(L):
    W = np.array([[I, sqrtlmdn,lmdn2],
                  [Z, I,  2*sqrtlmdn],
                  [Z, Z,     I      ]])
    Wfirst = np.array([[I, sqrtlmdn,lmdn2]])       # 左边界 MPO（只有第一行）

    Wlast = np.array([[lmdn2],[2*sqrtlmdn],[I]])        # 右边界 MPO（只有最后一列）

    MPO = [Wfirst] + ([W] * (L - 2)) + [Wlast]           # 组装完整 MPO 列表

    return MPO

#MPO加法
#---------------------------------------------------------------------------------------------------------------
def add_MPO(MPOA, MPOB):
    L = len(MPOA)
    assert len(MPOB) == L
    MPO_sum = []
    
    for i in range(L):
        WA = MPOA[i]
        WB = MPOB[i]
        d = WA.shape[2]          # 物理维数
        assert WA.shape[3] == d and WB.shape[2] == d and WB.shape[3] == d
        
        if i == 0:
            # 第一个张量：左维必须为1
            assert WA.shape[0] == 1 and WB.shape[0] == 1
            r1 = WA.shape[1]
            r2 = WB.shape[1]
            W_new = np.zeros((1, r1 + r2, d, d), dtype=WA.dtype)
            W_new[0, :r1, :, :] = WA[0, :, :, :]
            W_new[0, r1:, :, :] = WB[0, :, :, :]
            
        elif i == L - 1:
            # 最后一个张量：右维必须为1
            assert WA.shape[1] == 1 and WB.shape[1] == 1
            l1 = WA.shape[0]
            l2 = WB.shape[0]
            W_new = np.zeros((l1 + l2, 1, d, d), dtype=WA.dtype)
            W_new[:l1, 0, :, :] = WA[:, 0, :, :]
            W_new[l1:, 0, :, :] = WB[:, 0, :, :]
            
        else:
            # 中间张量：直接直和（左维和右维均相加），分块对角阵
            a, b = WA.shape[0], WA.shape[1]
            c, e = WB.shape[0], WB.shape[1]
            W_new = np.zeros((a + c, b + e, d, d), dtype=WA.dtype)
            W_new[:a, :b, :, :] = WA
            W_new[a:, b:, :, :] = WB
            
        MPO_sum.append(W_new)
        
    return MPO_sum

#构造粒子数惩罚项lmd*(N-N_target)^2，按完全平方写开
#---------------------------------------------------------------------------------------------------------------
def generate_MPO_penalty(L,excite):   #excite=0基态，+1粒子激发，-1空穴激发
    N_target=L+excite
    c=-2*lmd* N_target
    e= lmd*(N_target**2)
    W2 = np.array([[I, c*n],
                [Z, I]])
    W2first = np.array([[I,c*n]])       # 左边界 MPO（只有第一行）

    W2last = np.array([[c*n],[I]])        # 右边界 MPO（只有最后一列）

    MPO2 = [W2first] + ([W2] * (L - 2)) + [W2last]           # 组装完整 MPO 列表
    
    W3 = np.array([[I]]) 
    MPO3=[e*W3]+[W3]*(L-1)
    MPO1=generate_MPO_lmdN2(L)
    MPO2plus3=add_MPO(MPO2, MPO3)      
    MPO_penalty=add_MPO(MPO1,MPO2plus3)
    return MPO_penalty


#构造带惩罚的总哈密顿量
#---------------------------------------------------------------------------------------------------------------
def generate_MPO(t,L,excite):

    W = np.array([[I, ad, -t * a, h],         #Bose-hubbard模型原始哈密顿量MPO_H
                [Z, Z, Z, -t * a],            #形状(辅助维，辅助维数，物理维，物理维)=(4,4,d,d)
                [Z, Z, Z, ad],
                [Z, Z, Z,  I]])

    Wfirst = np.array([[I, ad, -t * a, h]])       # 左边界 MPO（只有第一行）

    Wlast = np.array([[h], [-t * a], [ad], [I]])        # 右边界 MPO（只有最后一列）

    MPO_H = [Wfirst] + ([W] * (L - 2)) + [Wlast]           # 组装完整 MPO 列表

    MPO_penalty=generate_MPO_penalty(L,excite)     

    MPO=add_MPO(MPO_H, MPO_penalty)          # 加上粒子数惩罚项
    return MPO,MPO_H                         # 除边界，MPO元素形状为(10,10,d,d),10=4+(3+2+1)



#真空
#---------------------------------------------------------------------------------------------------------------
def vaccum_E(W):
    E = np.zeros((W.shape[0], 1, 1))# 创建左边界真空环境
    E[0] = 1
    return E
def vaccum_F(W):
    F = np.zeros((W.shape[1], 1, 1))# 创建右边界真空环境 
    F[-1] = 1
    return F



#环境初始化
#---------------------------------------------------------------------------------------------------------------
def contract_from_left(W, A, E, B):
    Temp = np.einsum("sij,aik->sajk", A, E)
    Temp = np.einsum("sajk,abst->tbjk", Temp, W)
    return np.einsum("tbjk,tkl->bjl", Temp, B)

def contract_from_right(W, A, F, B):
    Temp = np.einsum("sij,bjl->sbil", A, F)
    Temp = np.einsum("sbil,abst->tail", Temp, W)
    return np.einsum("tail,tkl->aik", Temp, B)   

def initial_F(Alist, MPO, Blist):
    F = [vaccum_F(MPO[-1])]
    for i in range(len(MPO) - 1, 0, -1):
        F.append(contract_from_right(MPO[i], Alist[i], F[-1], Blist[i]))
    return F
    # 输出：列表 F，第一个元素是右真空（形状 (D_right_last, 1, 1)），后续每个元素是逐次收缩后的右环境

def initial_E(Alist, MPO, Blist):
    return [vaccum_E(MPO[0])]
    # 只返回左真空列表，因为从左边开始第一次dmrg扫描，初始左环境就是真空




#全收缩
#---------------------------------------------------------------------------------------------------------------
def expectation(AList, MPO, BList):
    E = [[[1]]]
    for i in range(0, len(MPO)):
        E = contract_from_left(MPO[i], AList[i], E, BList[i]) #递归
    return E[0][0][0]
    # 将整个 MPS 与 MPO 的张量网络完全收缩，得到哈密顿量的期望值 <psi| H |psi> 

# 粗粒化
#---------------------------------------------------------------------------------------------------------------
def coarse_grain_MPO(W, X):
    #将两个相邻的 MPO 张量 W 和 X 合并为一个粗粒化的 MPO 张量
    return np.reshape(np.einsum("abst,bcuv->acsutv", W, X),
                      [W.shape[0], X.shape[1], W.shape[2] * X.shape[2], W.shape[3] * X.shape[3]]) 

def coarse_grain_MPS(A, B):
    #将两个相邻的 MPS 张量 A 和 B 合并为一个粗粒化的波函数张量
    return np.reshape(np.einsum("sij,tjk->stik", A, B),
                      [A.shape[0] * B.shape[0], A.shape[1], B.shape[2]])




# 细粒化(svd)，截断(低秩近似)
#---------------------------------------------------------------------------------------------------------------
def fine_grain_MPS(A, dims):
    #将粗粒化的双格点波函数张量 A, 形状 (d*d, χL, χR)通过 SVD 分解回两个三阶张量：左格点左正则张量 U、奇异值数组 S、右格点右正则张量 V
    # 顺序：粗粒化 -> lancozs优化（只能作用于粗粒化的MPS）-> 细粒化，故并非多此一举（显得刚合起来又拆开）
    assert A.shape[0] == dims[0] * dims[1]
    Theta = np.transpose(np.reshape(A, dims + [A.shape[1], A.shape[2]]), (0, 2, 1, 3))

    # 先把要合并的索引transpose到一起，再用reshape合并，否则是一团乱麻！
    M = np.reshape(Theta, (dims[0] * A.shape[1], dims[1] * A.shape[2]))
    U, S, V = np.linalg.svd(M, full_matrices=0)
    U = np.reshape(U, (dims[0], A.shape[1], -1))
    V = np.transpose(np.reshape(V, (-1, dims[1], A.shape[2])), (1, 0, 2))

    return U, S, V

def truncate_MPS(U, S, V):
    
    total_norm2=np.sum(S**2)
    #print(f'奇异值平方和={total_norm2}')
    cum_norm2=np.cumsum(S**2)
    remain_number = np.searchsorted(cum_norm2, total_norm2 - target_trunc) + 1
    #searchsorted在有序数组 cumsum 中寻找第一个 >= total - target_trunc 的插入位置，返回该位置的索引，数量（维数）则+1
    
    m = min(max(remain_number, 1),D)  #自适应截断,1是为了健壮性，防止返回0维
    trunc = np.sum(S[m:]**2)/total_norm2
    #表示从索引 m 开始（包含）到末尾的所有元素，即实际保留态数 m
    S = S[0:m]
    U = U[:, :, 0:m]
    V = V[:, 0:m, :]
    return U, S, V, trunc, m




# 两点优化(lanczos)
#---------------------------------------------------------------------------------------------------------------
# scipy.sparse.linalg.eigsh 用lanczos方法求解稀疏矩阵本征值/矢，它不直接接受矩阵，而是接受一个 线性算子（linear operator） 对象，
# 该对象必须能够实现“矩阵乘以向量”的操作。只需要告诉它如何计算 H_eff|psi>而不需要把整个 H_eff矩阵 显式地构造出来。
class HamiltonianMultiply(sparse.linalg.LinearOperator):
    def __init__(self, E, W, F):
        # 存储左环境 E、合并后的 MPO W、右环境 F，并确定波函数的形状 req_shape 和总大小 size
        self.E = E
        self.W = W
        self.F = F
        self.dtype = np.dtype('d')
        self.req_shape = [W.shape[2], E.shape[1], F.shape[1]]
        self.size = self.req_shape[0] * self.req_shape[1] * self.req_shape[2]
        self.shape = [self.size, self.size]

    def _matvec(self, A):  #实现矩阵向量乘法 _matvec,是基类 `LinearOperator` 要求子类必须实现的方法,名称固定。
        # 输入向量 A 重塑为波函数张量，与环境张量收缩，得到 H_{eff}|A> $ 并展平返回。
        temp1 = np.einsum("aij,sik->ajsk", self.E, np.reshape(A, self.req_shape))
        temp2 = np.einsum("ajsk,abst->bjtk", temp1, self.W)
        R = np.einsum("bjtk,bkl->tjl", temp2, self.F)
        return np.reshape(R, -1) #上面收缩出的H_eff|A>仍为有三个索引的东西，得合并成一个才能传入eigsh
    

def optimize_two_sites(A, B, W1, W2, E, F, dir):
    
    # 粗粒化当前被优化格点对的 MPO 和 MPS 张量
    W = coarse_grain_MPO(W1, W2)
    AA = coarse_grain_MPS(A, B)
    
    # 实例化构造线性算子 H，即收缩出的有效H_eff
    H = HamiltonianMultiply(E, W, F)
    #Lanczos 求最小本征态
    Energy, V = sparse.linalg.eigsh(H, 1, v0=AA.ravel(), which='SA', maxiter=150, tol=1e-3)

    # 细粒化并截断
    AA = np.reshape(V[:, 0], H.req_shape)   #重塑回原来的三个索引的东西
    A, S, B = fine_grain_MPS(AA, [A.shape[0], B.shape[0]])
    A, S, B, trunc, m = truncate_MPS(A, S, B)
   
    # 根据扫描方向 'right' 或 'left' 将奇异值合并到右或左张量上，以移动正交中心
    if (dir == 'right'):
        B = np.einsum("ij,sjk->sik", np.diag(S), B)
    else:
        assert dir == 'left'
        A = np.einsum("sij,jk->sik", A, np.diag(S))

    # 返回能量、更新后的两个张量、截断误差和保留态数
    return Energy[0], A, B, trunc, m





# DMRG总流程
#---------------------------------------------------------------------------------------------------------------
def two_site_dmrg(MPS, MPO, sweeps):
    # MPS和MPO都是列表
    #初始化左/右环境
    E = initial_E(MPS, MPO, MPS) #左环境列表，初始只包含左真空
    F = initial_F(MPS, MPO, MPS) # 右环境列表，包含从右到左的所有环境
    F.pop()  # 移除 F 的最后一个元素（即包含所有格点的环境）
    #这样 F[-1] 就变成了包含格点 2..N-1 的环境，正好是优化第一对格点 (0,1) 所需的右环境（第一对格点左侧是空，右侧是格点 2..N-1）

    old_energy = expectation(MPS, MPO, MPS)   # 初始能量
    
    for sweep in range(0, int(sweeps / 2)):

        for i in range(0, len(MPS) - 2):
            # 从左到右优化所有相邻格点对，更新左环境
            Energy, MPS[i], MPS[i + 1], trunc, states = optimize_two_sites(MPS[i], MPS[i + 1],
                                                                            MPO[i], MPO[i + 1],
                                                                            E[-1], F[-1], 'right')
        
            E.append(contract_from_left(MPO[i], MPS[i], E[-1], MPS[i]))#更新左环境
            F.pop() 
            #移除右环境列表的最后一个元素，使 F[-1] 自动变成下一对格点所需的右环境（因为当前格点已被优化，应从右环境中去掉）。

        for i in range(len(MPS) - 2, 0, -1):# 左<-<-右,遍历MPS列表
            #e.g.N=5，(3,4)的情况在这里被优化了，防止把它连续优化两次
            Energy, MPS[i], MPS[i + 1], trunc, states = optimize_two_sites(MPS[i], MPS[i + 1],
                                                                            MPO[i], MPO[i + 1],
                                                                            E[-1], F[-1],'left')
           
            F.append(contract_from_right(MPO[i + 1], MPS[i + 1], F[-1], MPS[i + 1]))#更新右环境
            #关于为什么F还能加新东西。因为上面左->->右循环完后，F已经被pop空了，只剩下右真空
            E.pop()
            #同理，循环结束后，这里把E给pop空了，所以下一轮sweep才能append E

        new_energy = expectation(MPS, MPO, MPS)
        delta=abs(new_energy-old_energy)
        if delta<energy_tol:
            print(f"Sweep次数={sweep+1} 基态能量={new_energy:.12f} 变化={delta:.2e}")
            break
        old_energy = new_energy

    print(f"Sweep全部完成, 基态能量={new_energy:.12f} 变化={delta:.2e}")
    return MPS,new_energy




# 具体物理
#---------------------------------------------------------------------------------------------------------------
def compute_one_choice(t, L, excite): #计算一种(t,L,激发)的基态能量
    start = time.time()
    
    MPO ,MPO_H= generate_MPO(t, L,excite)
    N_MPO=generate_MPO_N(L)       #粒子数算符的MPO，用来监控粒子数是否守恒

    if excite == 0:
        MPS = generate_MPS(L)
    elif excite == 1:                #根据不同激发初始化MPS
        MPS = generate_MPS_p(L)
    else:
        MPS = generate_MPS_h(L)

    initial_N=expectation(MPS,N_MPO,MPS)
    print(f"*******************t={t:.2f}, L={L}, excite={excite}，初始总粒子数={initial_N}******************* ")
    
    Gs,E = two_site_dmrg(MPS, MPO, sweeps=12) #返回求出的基态MPS
    N = expectation(Gs, N_MPO, Gs)  #打印粒子数是否守恒

    E_H=expectation(Gs, MPO_H, Gs)    #E_ground=<Gs|H|Gs>
  
    waste_time = time.time() - start
    minutes = waste_time // 60
    seconds = waste_time % 60
    print(f"t={t:.2f}, L={L}, excite={excite} 基态求解完成！耗时: {int(minutes)}分{seconds:.1f}秒")
    print(f"基态总粒子数={N}")   
    print(f"真正基态能量，用H算的E=<H>={E_H}")
    print(f"------------------------------------------------------")
    return (t, L, excite, E_H)  #并行计算会打乱每个t, L, excite组合E的顺序！为了追踪回来，返回元组



def main():
    start = time.time()
    # 并行计算尤其是 multiprocessing 中，必须写，否则会出错。
    # 并行任务的创建和执行、主程序入口 放在它里面
    # 函数定义、类定义、常量声明等都可以放在外面。
    all_choices = [(t, L, excite) for L in L_list for excite in (0, 1, -1) for t in t_list]
    with mp.Pool(processes=8) as pool:
        # pool.starmap(func, iterable)将 iterable 中的每个元素解包后作为多个参数传给 func。
        results = pool.starmap(compute_one_choice, all_choices) 
        # results 是 pool.map 返回的一个列表，其中每个元素是 compute_one 函数的返回值。
        # 列表的顺序与 tasks 列表的顺序一致（因为用了 map）。

    #2026.4.16在此处报错：compute_one_choice返回元组，results的元素才能是元组，进行后续代码
    #情况1交互式环境（Jupyter Notebook、IPython）中运行
    #报错后，results 变量依然存在，可直接在交互式终端里修复后续代码，无需重新计算。

    #情况2：直接运行 python脚本.py，报错后进程会退出
    #内存中的数据全部丢失，results 无法恢复（除非提前写入了文件）。很遗憾只能重新运行，应该添加中间保存机制，避免未来再发生同样悲剧。
    
    import pickle
    # 保存 results 到文件,wb:write binary，以二进制写入模式打开文件, f:文件
    with open('energy_results_lobe2.pkl', 'wb') as f:
        pickle.dump(results, f)

    energy_dict = {}
    for t, L, excite, E in results:
        energy_dict[(t, L, excite)] = E  
        # 将 results 转换为字典，但元素无序！python不关心键—值对的存储顺序，而只跟踪键和值之间的关联关系。
    for i, t in enumerate(t_list):
        for j, L in enumerate(L_list):
            E0 = energy_dict[(t, L, 0)]
            Ep = energy_dict[(t, L, 1)]
            Eh= energy_dict[(t, L, -1)]
            
            mup_finite_list[i, j] = Ep - E0
            muh_finite_list[i, j] = E0-Eh


    #有限尺寸外推，固定t，把mu(t,L)对1/线性拟合，其截距（1/L=0，L=无穷）即热力学极限下的mu
    x = 1.0 / L_list 
    for i in range(len(t_list)):
        coeffs_p = np.polyfit(x, mup_finite_list[i,:], deg=1)# 截距=coeffs[1]  斜率=coeffs[0]
        mup_list[i]=coeffs_p[1]
        
        coeffs_h = np.polyfit(x, muh_finite_list[i,:], deg=1)
        muh_list[i]=coeffs_h[1]
    
        # 创建图形和坐标轴
    fig, ax = plt.subplots(figsize=(9, 6))

    # 绘制散点图（小方块）
    ax.scatter(t_list, mup_list, marker='s', s=36, edgecolors='black', linewidth=0.5,
            label='mu_p', alpha=0.8)
    ax.scatter(t_list, muh_list, marker='s', s=36, edgecolors='black', linewidth=0.5,
            label='mu_h', alpha=0.8)

    # 添加标签、标题、网格
    ax.set_xlabel('t', fontsize=12)
    ax.set_ylabel('mu', fontsize=12)
    ax.set_title('BH', fontsize=14)

    # ========== 添加自动刻度线 ==========
    from matplotlib.ticker import AutoLocator, AutoMinorLocator
    # 主刻度线：自动选择合适位置
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(AutoLocator())
    # 次刻度线：自动在相邻主刻度之间划分
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # 调整刻度线样式（长度、宽度、方向）
    ax.tick_params(which='major', length=10, width=1.5, direction='in')
    ax.tick_params(which='minor', length=6, width=1, direction='in')
    # =================================
    # 显示网格（普通网格线）
    ax.grid(True, linestyle='--', alpha=0.7)

    ax.legend()
    fig.tight_layout()

    # 保存高清图片
    plt.savefig('KT_transition.png', dpi=300)
    #py脚本下，show()运行的图并没被保存到任何地方，也没显示出来，而是在内存中被创建然后被丢弃了。
    
    waste_time = time.time() - start
    minutes = waste_time // 60
    seconds = waste_time % 60
    print(f"BH_DMRG全流程完成！总耗时: {int(minutes)}分{seconds:.1f}秒")


if __name__ == "__main__":
    main()
