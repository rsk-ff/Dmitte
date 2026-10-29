# %%
import numpy as np
from scipy.linalg import expm

from dmitte import para_constant
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant


# 迁移率矩阵的列顺序, 与 transfer_equotions() 的解包顺序一致
K_NAMES = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
           'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh',
           'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh',
           'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh', 'kbh_ba', 'kba_bh', 'kba_bo', 'kba_fo']
K_INDEX = {name: i for i, name in enumerate(K_NAMES)}
# 库室顺序, 与 transfer_equotions() 的 y 一致; ba (叶片新鲜同化物 OBT) 放在最后, 前 10 个的编号不变
A_NAMES = ['a1', 's0', 'a2', 's1', 's2', 's3', 'bh', 'bo', 'fh', 'fo', 'ba']
A_INDEX = {name: i for i, name in enumerate(A_NAMES)}


def leaf_exchange_rate(plantfx):
    '''
    叶片 -> 大气的 HTO 交换速率 kbh_a2, h-1

    与吸收项 ka2_bh = VDPF / ML 共用同一个总导度 g = VDPF = 1/(RAM+RB+RC1):
        吸收通量 = g * rho_v   * C_air            (rho_v: 空气水汽密度, C_air: 水汽比活度)
        释放通量 = g * rho_sat * C_leaf / ALF_HTO
    稳态时 C_leaf = ALF_HTO * RH * C_air。交换与净蒸腾无关, 夜间蒸腾趋于 0 时照常进行。
    rho_sat 取气温下的饱和水汽密度 ROW (g/m3), plant_w 为植物水量 (kg/m2)。
    '''
    return plantfx.VDPF * (plantfx.ROW / 1000) / (para_constant.ALF_HTO * plantfx.plant_w) * 3600


ROOT_LAYER_SHARE = (0.2, 0.4, 0.4)       # 根系吸水来自土壤第 1、2、3 层的份额


def net_transpiration(airfx, plantfx):
    '''
    净蒸腾 E = g * (rho_sat - rho_v), kg m-2 h-1, 与叶片-大气 HTO 交换共用导度 g = VDPF:
    E = 叶片水汽释放 (kbh_a2 * ALF_HTO * W_bh) - 叶片水汽吸收 (ka2_bh * W_air); 结露 (E < 0) 时取 0。
    '''
    gross_out = plantfx.VDPF * plantfx.ROW / 1000 * 3600
    gross_in = plantfx.VDPF / airfx.ML * 3600 * airfx.atm_h * 9
    return np.maximum(gross_out - gross_in, 0.0)


def _rate_of_change(w):
    '''逐时水量序列的变化率, kg m-2 h-1'''
    w = np.asarray(w, dtype=np.float64)
    return np.gradient(w) if w.size > 1 else np.zeros_like(w)


def plant_water_routes(airfx, soilfx, plantfx, kbh_so, kbh_fh, fruit_share):
    '''
    根系吸水和植物-果实液流, h-1, 使植物水 (bh) 和果实水 (fh) 的水量平衡闭合:
        根系总吸水 U = 净蒸腾 E + kbh_so * W_bh + d(W_bh + W_fh)/dt    (kg m-2 h-1, 不小于 0)
    U 按 ROOT_LAYER_SHARE 取自三层土壤, 其中 fruit_share 直接进入果实 (块茎), 其余进入植物水。
    植物 -> 果实液流取 kbh_fh * W_bh, 果实 -> 植物回流取余量, 使果实水的净流入等于 dW_fh/dt;
    余量为负时改为加大植物 -> 果实液流。fruit_share=None 表示没有果实库室 (叶菜)。

    Return: ks1_bh, ks2_bh, ks3_bh, ks1_fh, ks2_fh, ks3_fh, kbh_fh, kfh_bh
    '''
    n = len(airfx.U10)
    full = lambda v: np.broadcast_to(np.asarray(v, dtype=np.float64), (n,)).copy()
    has_fruit = fruit_share is not None
    w_bh = full(plantfx.plant_w)
    w_fh = full(plantfx.friut_w) if has_fruit else np.zeros(n)
    dw_bh, dw_fh = _rate_of_change(w_bh), _rate_of_change(w_fh)
    uptake = np.maximum(net_transpiration(airfx, plantfx) + kbh_so * w_bh + dw_bh + dw_fh, 0.0)
    to_fruit = (fruit_share if has_fruit else 0.0) * uptake
    f_bf = full(kbh_fh) * w_bh if has_fruit else np.zeros(n)
    f_fb = f_bf + to_fruit - dw_fh
    f_bf, f_fb = f_bf + np.maximum(-f_fb, 0), np.maximum(f_fb, 0)
    k_bf = np.divide(f_bf, w_bh, out=np.zeros(n), where=w_bh > 0)
    k_fb = np.divide(f_fb, w_fh, out=np.zeros(n), where=w_fh > 0)
    soil_w = (soilfx.soil1w, soilfx.soil2w, soilfx.soil3w)
    ks_bh = [share * (uptake - to_fruit) / w for share, w in zip(ROOT_LAYER_SHARE, soil_w)]
    ks_fh = [share * to_fruit / w for share, w in zip(ROOT_LAYER_SHARE, soil_w)]
    return (*ks_bh, *ks_fh, k_bf, k_fb)


def _split_darcy(k_down, k_up, w_upper, w_lower):
    '''
    把可能为负的向下迁移率 k_down (上层 -> 下层) 拆成两个非负速率。
    k_down < 0 表示水向上流: 通量 |k_down| * w_upper 改由下层提供, 加到向上的速率 k_up 上。
    '''
    k_down = np.asarray(k_down, dtype=np.float64)
    k_up = k_up + np.maximum(-k_down, 0) * w_upper / w_lower
    return np.maximum(k_down, 0), k_up


def _assemble(k_list, n_t_step, plantfx):
    '''
    把各迁移率组装成 (n_t_step, len(K_NAMES)) 的迁移率矩阵, 并处理还不存在的植物库室:
    - 植物水量为 0 (出苗前) 的时刻, 进出植物和果实库室的迁移率都为 0;
    - 果实/块茎水量为 0 (形成前) 的时刻, 根系吸水全部进入植物水, 不向果实迁移。
    '''
    k = np.array([np.full(n_t_step, v, dtype=np.float64) if np.isscalar(v) else np.asarray(v, dtype=np.float64)
                  for v in k_list]).T
    with np.errstate(invalid='ignore'):
        no_plant = ~(np.broadcast_to(plantfx.plant_w, (n_t_step,)) > 0)
        no_fruit = ~(np.broadcast_to(plantfx.friut_w, (n_t_step,)) > 0)
    for layer in '123':
        to_bh, to_fh = K_INDEX['ks%s_bh' % layer], K_INDEX['ks%s_fh' % layer]
        k[no_fruit, to_bh] += k[no_fruit, to_fh]
    fruit_rates = [K_INDEX[n] for n in ('ks1_fh', 'ks2_fh', 'ks3_fh', 'kbh_fh', 'kfh_bh', 'kbh_fo', 'kba_fo')]
    plant_rates = fruit_rates + [K_INDEX[n] for n in ('ks1_bh', 'ks2_bh', 'ks3_bh', 'ka2_bh', 'kbh_a2',
                                                       'kbh_so', 'kbh_bo', 'kbo_bh', 'kbh_ba', 'kba_bh', 'kba_bo')]
    k[np.ix_(no_fruit, fruit_rates)] = 0
    k[np.ix_(no_plant, plant_rates)] = 0
    return k


def check_rates(k_array):
    '''迁移率必须有限且非负, 否则抛出 ValueError 并列出有问题的迁移率。'''
    k = np.asarray(k_array, dtype=np.float64)
    names = K_NAMES if k.shape[-1] == len(K_NAMES) else [str(j) for j in range(k.shape[-1])]
    problems = []
    for j, name in enumerate(names):
        n_bad = int((~np.isfinite(k[..., j])).sum())
        n_neg = int((k[..., j] < 0).sum())
        if n_bad:
            problems.append('%s: %d non-finite' % (name, n_bad))
        if n_neg:
            problems.append('%s: %d negative' % (name, n_neg))
    if problems:
        raise ValueError('Transfer rates must be finite and >= 0 (' + '; '.join(problems) + ')')


def water_balance_residual(k_array, airfx, soilfx, plantfx, alpha=para_constant.ALF_HTO):
    '''
    均匀比活度检验 (水量闭合检验), h-1

    令所有水库室的比活度相同 (A_i = W_i)、不计衰变和 OBT 生成/分解、去掉叶片释放中的同位素分馏
    (kbh_a2 乘以 alpha; kbh_a2 不含分馏因子的变体传 alpha=1), 此时 dA_i/dt / A_i 就是各迁移率隐含的
    净水通量 / 水量。水量平衡闭合时它应等于水量的相对变化率 dlnW/dt, 通常远小于库室的周转速率 k_out;
    偏差就是均匀比活度下出现的虚假比活度变化。大气作为边界条件, 只计它流向土壤和植物的通量。

    Return
    ---
    {库室: {'residual', 'k_out', 'dlnW_dt'}}, 库室为 s1, s2, s3, bh, fh, 各项为逐时数组
    '''
    k = np.array(k_array, dtype=np.float64, copy=True)
    n = len(k)
    for name in ('kbh_bo', 'kbh_fo', 'kbo_bh', 'kbh_ba', 'kba_bh', 'kba_bo', 'kba_fo'):
        k[:, K_INDEX[name]] = 0
    k[:, K_INDEX['kbh_a2']] *= alpha
    W = {'a2': airfx.atm_h * 9, 's1': soilfx.soil1w, 's2': soilfx.soil2w, 's3': soilfx.soil3w,
         'bh': plantfx.plant_w, 'fh': plantfx.friut_w}
    W = {c: np.broadcast_to(np.asarray(w, dtype=np.float64), (n,)) for c, w in W.items()}
    zero = np.zeros(n)
    y = [zero, zero, W['a2'], W['s1'], W['s2'], W['s3'], W['bh'], zero, W['fh'], zero, zero]
    dy = transfer_equotions(0, y, k.T, 0.0)
    position = {'s1': 3, 's2': 4, 's3': 5, 'bh': 6, 'fh': 8}
    outflows = {'s1': ('ks1_a2', 'ks1_bh', 'ks1_fh', 'ks1_s2'),
                's2': ('ks2_s1', 'ks2_bh', 'ks2_fh', 'ks2_s3'),
                's3': ('ks3_s3', 'ks3_bh', 'ks3_fh', 'ks3_s2'),
                'bh': ('kbh_a2', 'kbh_so', 'kbh_fh'),
                'fh': ('kfh_bh',)}
    out = {}
    with np.errstate(divide='ignore', invalid='ignore'):
        for c, i in position.items():
            out[c] = {'residual': dy[i] / W[c],
                      'k_out': sum(k[:, K_INDEX[name]] for name in outflows[c]),
                      'dlnW_dt': np.gradient(np.log(W[c])) if n > 1 else zero}
    return out


def transfer_rates(char: str,
                   plant_type: str,
                   airfx: Calc_air,
                   soilfx: Calc_soil,
                   plantfx: Calc_plant,
                   percolation: str = para_constant.PERCOLATION,
                   kbh_fh: float = para_constant.KBH_FH) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG'
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    percolation: 土壤层间下渗, 'constant' 为常数通量 (第 1->2 层 0.181*9/24, 第 2->3 层 0.072*9/24 kg m-2 h-1),
                 'darcy' 为 calc_soil_ab 的达西通量 (向上流时由 _split_darcy 拆成下层 -> 上层)
    kbh_fh: 植物水 -> 果实水的交换速率, h-1; 生长所需的净液流由 plant_water_routes 另加

    Return:
    ---
    k_array: 迁移率矩阵, (小时数, len(K_NAMES)), 列顺序见 K_NAMES
    '''
    ka2_s1 = (soilfx.VDSO / airfx.ML * 3600) + (airfx.rainfall / 9) / airfx.atm_h
    ks3_s2 = ks3_s3 = 3.42E-4
    ks2_s1 = ks3_s3 * soilfx.soil3w / soilfx.soil2w
    if percolation == 'darcy':
        ks1_s2, ks2_s1 = _split_darcy(soilfx.Va_b1 / 100, ks2_s1, soilfx.soil1w, soilfx.soil2w)
        ks2_s3, ks3_s2 = _split_darcy(soilfx.Va_b2 / 150, ks3_s2, soilfx.soil2w, soilfx.soil3w)
    elif percolation == 'constant':
        ks1_s2 = 0.181 * 9 / soilfx.soil1w / 24
        ks2_s3 = 0.072 * 9 / soilfx.soil2w / 24
    else:
        raise ValueError("Unsupported percolation: '{}'.".format(percolation))
    kbh_so = 3E-4
    ka2_bh = plantfx.VDPF / airfx.ML * 3600          # 降雨不截留在植物上, 只进入土壤 (ka2_s1)
    kbh_a2 = leaf_exchange_rate(plantfx)             # 叶片-大气 HTO 对称交换
    kbo_bh = para_constant.KBO_BH

    if char.upper() == 'HT':
        ka1_a1 = 0.693
        ka2_a2 = 0
        ka1_a2 = 4E-3
        ka1_s0 = 4.5e-4 / airfx.ML * 3600
        ks0_s1 = 0.2
        ks0_s2 = 0.02
        ks0_s3 = 0.02
        ks0_a1 = ks1_a2 = soilfx.ESOIL / soilfx.soil1w * (0.5+0.1)

    elif char.upper() == 'HTO':
        ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
        ka2_a2 = 0.693
        ks1_a2 = soilfx.ESOIL / soilfx.soil1w * (1+0.1)

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))

    # OBT: 叶片水 -> 同化物库 ba -> 茎叶 bo / 果实 fo / 回到叶片水 (calc_transOBT); kbh_bo 只补茎叶 OBT 的周转。
    # 果实 OBT 全部经同化物库生成, 不再从叶片水直接生成 (kbh_fo = 0)
    kbh_bo = plantfx.TROBT
    kbh_fo = 0
    kbh_ba, kba_bh, kba_bo, kba_fo = plantfx.TR_ASSIM, plantfx.KBA_BH, plantfx.KBA_BO, plantfx.KBA_FO
    # fruit_share: 根系吸水中直接进入果实 (块茎) 的份额; None 表示没有果实库室
    if plant_type.upper() == 'LEAFY_PLANT':
        fruit_share = None
    elif plant_type.upper() == 'CEREAL':
        fruit_share = 0
    elif plant_type.upper() == 'ROOT_VEG':
        fruit_share = 1 / 6
    else:
        raise ValueError("Unsupported plant type: '{}'.".format(plant_type))

    ks1_bh, ks2_bh, ks3_bh, ks1_fh, ks2_fh, ks3_fh, kbh_fh, kfh_bh = plant_water_routes(
        airfx, soilfx, plantfx, kbh_so, kbh_fh, fruit_share)

    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2,
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh,
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh,
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh, kbh_ba, kba_bh, kba_bo, kba_fo]

    return _assemble(k_list, len(airfx.U10), plantfx)


def transfer_equotions(t, y, transfer_rates, LAMBDA_T: float):
    '''
    氚转移微分方程组
    
    Parameters:
    ---
    y, t: 形式需要, y 为待求解的变量向量, t 表示时间
    transfer_rates: d-1 or h-1, 迁移率列表或向量
    LAMBDA_T: d-1 or h-1, 氚的衰变常数, 单位与迁移率保持一致
    '''
    A_a1, A_s0, A_a2, A_s1, A_s2, A_s3, A_bh, A_bo, A_fh, A_fo, A_ba = y
    ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, ka2_s1, \
        ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, ks3_bh, ka2_a2, \
        kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, kbh_fh, kbh_fo, ks1_fh, \
        ks2_fh, ks3_fh, kbh_ba, kba_bh, kba_bo, kba_fo = transfer_rates

    
    dA_a1_dt = ks0_a1 * A_s0 - (ka1_a1 + ka1_a2 + ka1_s0 + LAMBDA_T) * A_a1
    dA_s0_dt = ka1_s0 * A_a1 - (ks0_a1 + ks0_s1 + ks0_s2 + ks0_s3 + LAMBDA_T) * A_s0
    dA_a2_dt = ka1_a2 * A_a1 + ks1_a2 * A_s1 + kbh_a2 * A_bh - (ka2_a2 + ka2_bh + ka2_s1 + LAMBDA_T) * A_a2
    dA_s1_dt = ks0_s1 * A_s0 + ka2_s1 * A_a2 + ks2_s1 * A_s2 + kbh_so * A_bh - (ks1_a2 + ks1_bh + ks1_fh + ks1_s2 + LAMBDA_T) * A_s1
    dA_s2_dt = ks0_s2 * A_s0 + ks1_s2 * A_s1 + ks3_s2 * A_s3 - (ks2_s1 + ks2_bh + ks2_fh + ks2_s3 + LAMBDA_T) * A_s2
    dA_s3_dt = ks0_s3 * A_s0 + ks2_s3 * A_s2 - (ks3_s3 + ks3_bh + ks3_fh + ks3_s2 + LAMBDA_T) * A_s3
    dA_bh_dt = ka2_bh * A_a2 + ks1_bh * A_s1 + ks2_bh * A_s2 + ks3_bh * A_s3 + kbo_bh * A_bo + kfh_bh * A_fh + kba_bh * A_ba - (kbh_a2 + kbh_so + kbh_bo + kbh_fh + kbh_fo + kbh_ba + LAMBDA_T) * A_bh
    dA_bo_dt = kbh_bo * A_bh + kba_bo * A_ba - (kbo_bh + LAMBDA_T) * A_bo
    dA_fh_dt = ks1_fh * A_s1 + ks2_fh * A_s2 + ks3_fh * A_s3 + kbh_fh * A_bh - (kfh_bh + LAMBDA_T) * A_fh
    dA_fo_dt = kbh_fo * A_bh + kba_fo * A_ba - LAMBDA_T * A_fo
    dA_ba_dt = kbh_ba * A_bh - (kba_bh + kba_bo + kba_fo + LAMBDA_T) * A_ba

    return dA_a1_dt, dA_s0_dt, dA_a2_dt, dA_s1_dt, dA_s2_dt, dA_s3_dt, dA_bh_dt, dA_bo_dt, dA_fh_dt, dA_fo_dt, dA_ba_dt


def solve_con(char, init_con, k_array, LAMBDA_T):
    '''
    计算各库室浓度

    Parameters:
    ---
    char: 'HT' or 'HTO'
    init_con: 初始浓度, Bq/m2
    k_array: 迁移率矩阵
    LAMBDA_T: 氚的衰变常数

    Return:
    ---
    As: 库室浓度矩阵, (小时数 + 1, len(A_NAMES)), 列顺序见 A_NAMES
    '''
    n_comp = len(A_NAMES)
    eye = np.eye(n_comp)

    def next_con(Ai, ki):
        '''
        定义单步求解过程, 请勿外部调用。
        方程组对 A 线性, 一小时内迁移率不变, 所以一小时的解就是生成矩阵 (逐列取自 transfer_equotions) 的矩阵指数;
        果实形成初期等时刻迁移率可达 1e3-1e4 h-1, 显式 RK45 在这些时刻误差可达 100 %。
        '''
        generator = np.column_stack([transfer_equotions(0, eye[j], ki, LAMBDA_T) for j in range(n_comp)])
        return expm(generator) @ Ai
    # Preserve caller-owned transfer rates; the first-step adjustments are solver-local.
    k_array = np.array(k_array, dtype=np.float64, copy=True)
    check_rates(k_array)
    char = char.upper()
    k_array[0, 0] = 0
    k_array[0, 17] = 0

    n_t_step = len(k_array)
    As = np.zeros((n_t_step + 1, n_comp))

    if char == 'HT':
        As[0, A_INDEX['a1']] = init_con

        As[1] = next_con(As[0], k_array[0])

        As[1, 0] = 0
        As[1, 2] = 0
        
        for i in range(1, n_t_step):
            As[i + 1] = next_con(As[i], k_array[i])

    elif char == 'HTO':
        As[0, A_INDEX['a2']] = init_con

        As[1] = next_con(As[0], k_array[0])

        As[1, 0] = 0
        As[1, 2] = 0
        
        for i in range(1, n_t_step):
            As[i + 1] = next_con(As[i], k_array[i])

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))

    return As
