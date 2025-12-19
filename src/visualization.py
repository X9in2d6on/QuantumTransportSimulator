import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 加载结果数据
def load_data(filename):
    """加载模拟器生成的CSV结果文件"""
    print(f"加载数据文件: {filename}")
    df = pd.read_csv(filename)
    return df

# 绘制电势分布 2D 截面图
def plot_potential_2d(df, z_slice=5.0):
    """绘制特定z位置的电势分布2D截面图"""
    print(f"绘制z={z_slice}nm平面的电势分布")
    
    # 查找最接近指定z值的切片
    df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
    
    if len(df_slice) == 0:
        print(f"警告: 没有找到z={z_slice}nm附近的数据，尝试使用可用的z值")
        z_values = sorted(df['z(nm)'].unique())
        if len(z_values) > 0:
            z_slice = z_values[len(z_values)//2]  # 使用中间的z值
            print(f"改用z={z_slice}nm")
            df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
        else:
            print("错误: 没有可用的z值数据")
            return
    
    # 提取唯一的x和y坐标
    x_values = np.sort(df_slice['x(nm)'].unique())
    y_values = np.sort(df_slice['y(nm)'].unique())
    
    if len(x_values) == 0 or len(y_values) == 0:
        print("错误: 没有足够的x或y值数据")
        return
    
    # 创建网格
    X, Y = np.meshgrid(x_values, y_values)
    potential = np.zeros_like(X)
    
    # 填充电势值
    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            row = df_slice[(df_slice['x(nm)'] == x) & (df_slice['y(nm)'] == y)]
            if not row.empty:
                potential[j, i] = row['potential(V)'].values[0]
    
    # 绘制电势分布
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(X, Y, potential, 50, cmap='viridis')
    plt.colorbar(contour, label='Potential (V)')
    plt.xlabel('X Position (nm)')
    plt.ylabel('Y Position (nm)')
    plt.title(f'Potential Distribution at Z = {z_slice} nm')
    plt.savefig(f'potential_2d_z{z_slice}.png', dpi=300)
    plt.close()

# 绘制沿通道方向的电势分布
def plot_potential_along_channel(df, y_middle=None, z_middle=None):
    """绘制沿着通道方向的电势分布"""
    print("绘制沿通道方向的电势分布")
    
    # 如果未指定y和z，使用中间值
    if y_middle is None:
        y_middle = (df['y(nm)'].max() + df['y(nm)'].min()) / 2
    if z_middle is None:
        z_middle = (df['z(nm)'].max() + df['z(nm)'].min()) / 2
    
    # 查找最接近指定y和z值的切片
    df_slice = df[(np.abs(df['y(nm)'] - y_middle) < 0.5) & (np.abs(df['z(nm)'] - z_middle) < 0.5)]
    
    if len(df_slice) == 0:
        print("警告: 没有找到通道数据，尝试使用所有数据的中间位置")
        y_middle = df['y(nm)'].median()
        z_middle = df['z(nm)'].median()
        df_slice = df[(np.abs(df['y(nm)'] - y_middle) < 1.0) & (np.abs(df['z(nm)'] - z_middle) < 1.0)]
    
    if len(df_slice) == 0:
        print("错误: 仍然没有找到通道数据")
        return
    
    # 按x排序
    df_slice = df_slice.sort_values('x(nm)')
    
    # 绘制电势分布
    plt.figure(figsize=(12, 6))
    plt.plot(df_slice['x(nm)'], df_slice['potential(V)'], 'b-', linewidth=2)
    plt.xlabel('X Position (nm)')
    plt.ylabel('Potential (V)')
    plt.title(f'Potential Along Channel Center (y={y_middle:.1f}nm, z={z_middle:.1f}nm)')
    plt.grid(True)
    plt.savefig('potential_along_channel.png', dpi=300)
    plt.close()

# 绘制电子密度分布
def plot_electron_density(df, z_slice=5.0):
    """绘制特定z位置的电子密度分布"""
    print(f"绘制z={z_slice}nm平面的电子密度分布")
    
    # 查找最接近指定z值的切片
    df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
    
    if len(df_slice) == 0:
        print(f"警告: 没有找到z={z_slice}nm附近的数据，尝试使用可用的z值")
        z_values = sorted(df['z(nm)'].unique())
        if len(z_values) > 0:
            z_slice = z_values[len(z_values)//2]
            print(f"改用z={z_slice}nm")
            df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
        else:
            print("错误: 没有可用的z值数据")
            return
    
    # 提取唯一的x和y坐标
    x_values = np.sort(df_slice['x(nm)'].unique())
    y_values = np.sort(df_slice['y(nm)'].unique())
    
    if len(x_values) == 0 or len(y_values) == 0:
        print("错误: 没有足够的x或y值数据")
        return
    
    # 创建网格
    X, Y = np.meshgrid(x_values, y_values)
    density = np.zeros_like(X)
    
    # 填充电子密度值 (使用对数比例)
    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            row = df_slice[(df_slice['x(nm)'] == x) & (df_slice['y(nm)'] == y)]
            if not row.empty:
                # 取对数，避免数值差异过大
                density[j, i] = np.log10(max(row['electron_density(cm-3)'].values[0], 1.0))
    
    # 绘制电子密度分布
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(X, Y, density, 50, cmap='plasma')
    plt.colorbar(contour, label='Electron Density (log₁₀ cm⁻³)')
    plt.xlabel('X Position (nm)')
    plt.ylabel('Y Position (nm)')
    plt.title(f'Electron Density Distribution at Z = {z_slice} nm')
    plt.savefig(f'electron_density_z{z_slice}.png', dpi=300)
    plt.close()

# 绘制量子修正因子分布
def plot_quantum_factor(df, z_slice=5.0):
    """绘制特定z位置的量子修正因子分布"""
    print(f"绘制z={z_slice}nm平面的量子修正因子分布")
    
    # 查找最接近指定z值的切片
    df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
    
    if len(df_slice) == 0:
        print(f"警告: 没有找到z={z_slice}nm附近的数据，尝试使用可用的z值")
        z_values = sorted(df['z(nm)'].unique())
        if len(z_values) > 0:
            z_slice = z_values[len(z_values)//2]
            print(f"改用z={z_slice}nm")
            df_slice = df[np.abs(df['z(nm)'] - z_slice) < 0.1]
        else:
            print("错误: 没有可用的z值数据")
            return
    
    # 提取唯一的x和y坐标
    x_values = np.sort(df_slice['x(nm)'].unique())
    y_values = np.sort(df_slice['y(nm)'].unique())
    
    if len(x_values) == 0 or len(y_values) == 0:
        print("错误: 没有足够的x或y值数据")
        return
    
    # 创建网格
    X, Y = np.meshgrid(x_values, y_values)
    qfactor = np.zeros_like(X)
    
    # 填充量子修正因子值
    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            row = df_slice[(df_slice['x(nm)'] == x) & (df_slice['y(nm)'] == y)]
            if not row.empty:
                qfactor[j, i] = row['quantum_factor'].values[0]
    
    # 绘制量子修正因子分布
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(X, Y, qfactor, 50, cmap='coolwarm')
    plt.colorbar(contour, label='Quantum Correction Factor')
    plt.xlabel('X Position (nm)')
    plt.ylabel('Y Position (nm)')
    plt.title(f'Quantum Factor Distribution at Z = {z_slice} nm')
    plt.savefig(f'quantum_factor_z{z_slice}.png', dpi=300)
    plt.close()

# 绘制I-V特性曲线
def plot_iv_curves(filename='iv_curves.csv'):
    """绘制I-V特性曲线"""
    print(f"绘制I-V特性曲线，数据源：{filename}")
    
    # 加载I-V数据
    iv_data = pd.read_csv(filename)
    
    # 准备绘图
    plt.figure(figsize=(12, 8))
    
    # 获取所有不同的栅极电压
    vgs_values = sorted(iv_data['Vgs(V)'].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(vgs_values)))
    
    # 针对每个Vgs绘制Id-Vds曲线
    for i, vgs in enumerate(vgs_values):
        vgs_data = iv_data[iv_data['Vgs(V)'] == vgs]
        vgs_data = vgs_data.sort_values('Vds(V)')
        
        # 将电流转换为微安
        current_uA = vgs_data['Ids(A)'] * 1e6
        
        plt.plot(vgs_data['Vds(V)'], current_uA, 'o-', color=colors[i], linewidth=2, 
                 label=f'Vgs = {vgs} V')
    
    plt.xlabel('Drain-Source Voltage Vds (V)')
    plt.ylabel('Drain current Ids (μA)')
    plt.title('Nano FinFET Output Characterization Curves')
    plt.legend()
    plt.grid(True)
    plt.savefig('output_characteristics.png', dpi=300)
    plt.close()
    
    # 绘制转移特性 (Id-Vgs)
    plt.figure(figsize=(12, 8))
    
    # 获取所有不同的漏极电压
    vds_values = sorted(iv_data['Vds(V)'].unique())
    colors = plt.cm.plasma(np.linspace(0, 1, len(vds_values)))
    
    # 针对每个Vds绘制Id-Vgs曲线
    for i, vds in enumerate(vds_values):
        vds_data = iv_data[iv_data['Vds(V)'] == vds]
        vds_data = vds_data.sort_values('Vgs(V)')
        
        # 将电流转换为微安
        current_uA = vds_data['Ids(A)'] * 1e6
        
        plt.plot(vds_data['Vgs(V)'], current_uA, 'o-', color=colors[i], linewidth=2,
                 label=f'Vds = {vds:.1f} V')
    
    plt.xlabel('Gate Source Voltage Vgs (V)')
    plt.ylabel('Drain current Ids (μA)')
    plt.title(' Nanometer FinFET Transfer Characterization Curves')
    plt.legend()
    plt.grid(True)
    plt.savefig('transfer_characteristics.png', dpi=300)
    
    # 绘制对数比例的转移特性曲线
    plt.figure(figsize=(12, 8))
    
    # 针对每个Vds绘制log(Id)-Vgs曲线
    for i, vds in enumerate(vds_values):
        vds_data = iv_data[iv_data['Vds(V)'] == vds]
        vds_data = vds_data.sort_values('Vgs(V)')
        
        # 将电流转换为微安并取对数
        current_uA = vds_data['Ids(A)'] * 1e6
        current_uA = current_uA.clip(lower=1e-3)  # 避免负值或零值
        
        plt.semilogy(vds_data['Vgs(V)'], current_uA, 'o-', color=colors[i], linewidth=2,
                     label=f'Vds = {vds:.1f} V')
    
    plt.xlabel('Gate Source Voltage Vgs (V)')
    plt.ylabel('Drain current Ids (μA) [Log Scale]')
    plt.title('Nano FinFET Logarithmic Transfer Characterization Curves')
    plt.legend()
    plt.grid(True)
    plt.savefig('log_transfer_characteristics.png', dpi=300)
    plt.close()

# 主函数
def main():
    # 结果数据文件
    results_file = "nano_mosfet_results.csv"
    iv_file = "iv_curves.csv"
    
    # 加载数据
    try:
        df = load_data(results_file)
        print(f"成功加载数据，共{len(df)}个数据点")
        
        # 绘制2D截面图
        plot_potential_2d(df, z_slice=5.0)
        
        # 绘制沿通道方向的电势
        plot_potential_along_channel(df)
        
        # 绘制电子密度分布
        plot_electron_density(df, z_slice=5.0)
        
        # 绘制量子修正因子分布
        plot_quantum_factor(df, z_slice=5.0)
        
        print("2D分布图绘制完成")
    except Exception as e:
        print(f"加载或绘制设备数据时出错: {e}")
    
    # 绘制I-V特性曲线
    try:
        plot_iv_curves(iv_file)
        print("I-V特性曲线绘制完成")
    except Exception as e:
        print(f"绘制I-V特性曲线时出错: {e}")
    
    print("所有可视化任务完成！")

if __name__ == "__main__":
    main()