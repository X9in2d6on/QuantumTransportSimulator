# QuantumTransportSimulator
纳米尺度FinFET器件三维量子输运模拟框架使用手册
 1. 简介
本框架实现了一个基于纳米尺度FinFET器件的三维量子输运模拟系统，能够模拟考虑量子效应的电子输运过程。框架采用C++编写，具有高性能和可扩展性，支持以下功能：
 三维设备结构的精确模拟
 量子效应修正的自适应调整
 自洽泊松薛定谔方程求解
 考虑界面态影响的载流子动力学
 漂移扩散模型与量子效应的耦合
 器件特性曲线的生成与分析

 2. 理论基础
 2.1 量子修正因子模型
为了在传统的漂移扩散模型中考虑量子效应，本框架采用量子修正因子方法。其核心公式为：
```
γ = γ₀ * exp(α * |∇φ|)
```

其中：
 γ 是量子修正因子，用于调整载流子密度计算
 γ₀ 是基准修正系数（默认值3.0）
 α 是与电场强度相关的参数（默认值6e9）
 |∇φ| 是局部电势梯度的绝对值

 2.2 载流子密度计算
考虑量子效应后，载流子密度计算公式修正为：

```
n = ni * exp(q*φ/(kT*γ))
p = ni * exp(q*φ/(kT*γ))
```
通过引入修正因子γ，能够模拟量子限制效应导致的载流子分布变化。

 2.3 泊松方程
系统通过求解以下泊松方程计算电势分布：
```
∇²φ = q/ε * (p  n + Nd  Na + ρit)
```
其中ρit表示界面态电荷密度。

 2.4 电流计算
电流密度计算考虑了漂移和扩散两部分：
```
J_n = q*μ_n*n*E + q*D_n*∇n
J_p = q*μ_p*p*E  q*D_p*∇p
```
总电流通过在漏极边界上积分电流密度计算得到。

 3. 主要类与函数
 3.1 QuantumTransportSimulator 类
这是模拟框架的核心类，提供以下主要功能：
 构造函数
```cpp
QuantumTransportSimulator(
    double length,           // 设备长度 (m)
    double width,            // 设备宽度 (m)
    double height,           // 设备高度 (m)
    int nx, int ny, int nz,  // 各方向网格点数
    double effective_mass,   // 归一化有效质量 (m0单位)
    double oxide_thickness,  // 氧化层厚度 (m)
    double dielectric_constant // 相对介电常数
);
```

 电极设置
```cpp
// 设置源极和漏极
void setSourceDrain(
    double source_x_start, double source_x_end,
    double drain_x_start, double drain_x_end,
    double source_voltage, double drain_voltage
);

// 设置栅极
void setGate(
    double gate_x_start, double gate_x_end,
    double gate_y_start, double gate_y_end,
    double gate_z_position,
    double gate_voltage
);

// 设置衬底
void setSubstrate(
    double substrate_z_position, 
    double substrate_voltage
);
```
 掺杂设置

```cpp
// 设置细粒度掺杂分布
void setDoping(
    const Eigen::VectorXd& donor_concentration,
    const Eigen::VectorXd& acceptor_concentration
);

// 设置区域均匀掺杂
void setUniformDoping(
    double x_start, double x_end,
    double y_start, double y_end,
    double z_start, double z_end,
    double donor_conc, double acceptor_conc
);
```

 模拟与结果处理

```cpp
// 执行模拟
void simulate();

// 保存结果到文件
void saveResults(const std::string& filename);

// 获取漏电流
double getDrainCurrent() const;

// 获取电势分布
Eigen::VectorXd getPotential() const;

// 获取电子密度分布
Eigen::VectorXd getElectronDensity() const;
```

 模拟参数设置

```cpp
// 设置最大迭代次数
void setMaxIterations(int max_iter);

// 设置收敛容差
void setConvergenceTolerance(double tol);

// 设置量子修正因子参数
void setQuantumFactorParameters(double gamma0, double alpha);
```

 4. 使用示例
 4.1 基本器件模拟
以下是一个简单的FinFET器件模拟示例：

```cpp
include "QuantumTransportSimulator.h"

int main() {
    // 创建模拟器实例（所有尺寸单位：米）
    QuantumTransportSimulator simulator(
        35e9,   // 长度：35nm
        20e9,   // 宽度：20nm
        10e9,   // 高度：10nm
        70, 40, 30, // 网格分辨率
        0.25,    // 有效质量 (m0单位)
        1.2e9,  // 氧化层厚度：1.2nm
        3.9      // 二氧化硅介电常数
    );
    
    // 设置电极
    simulator.setSourceDrain(
        0, 5e9,           // 源极位置
        30e9, 35e9,      // 漏极位置
        0.0, 0.6           // 源极和漏极电压
    );
    
    simulator.setGate(
        5e9, 30e9,       // 栅极x范围
        0, 20e9,          // 栅极y范围
        10e9,             // 栅极z位置
        1.2                // 栅极电压
    );
    
    simulator.setSubstrate(0, 0.0);  // 衬底位置和电压
    
    // 设置掺杂（通道：P型，源漏：N型）
    simulator.setUniformDoping(
        5e9, 30e9, 0, 20e9, 0, 5e9,  // 通道区域
        0.0, 1e17                         // P型掺杂
    );
    
    simulator.setUniformDoping(
        0, 5e9, 0, 20e9, 0, 10e9,     // 源极区域
        5e19, 0.0                         // N型掺杂
    );
    
    simulator.setUniformDoping(
        30e9, 35e9, 0, 20e9, 0, 10e9, // 漏极区域
        5e19, 0.0                          // N型掺杂
    );
    
    // 执行模拟
    simulator.simulate();
    
    // 保存结果
    simulator.saveResults("results.csv");
    
    // 获取漏电流
    double current = simulator.getDrainCurrent();
    std::cout << "漏电流: " << current << " A" << std::endl;
    
    return 0;
}
```

 4.2 输出特性曲线生成
为了生成输出特性曲线，可以在不同的漏极电压和栅极电压下重复模拟，例如：
```cpp
// 栅极电压扫描
std::vector<double> gate_voltages = {0.4, 0.6, 0.8, 1.0, 1.2};
// 漏极电压扫描
std::vector<double> drain_voltages = {0.0, 0.2, 0.4, 0.6, 0.8};

// 结果存储
std::ofstream output("iv_curves.csv");
output << "Vgs(V),Vds(V),Ids(A)" << std::endl;

// 双重循环扫描
for (double vgs : gate_voltages) {
    for (double vds : drain_voltages) {
        // 创建模拟器并设置参数...
        
        // 设置电压
        simulator.setSourceDrain(0, 5e9, 30e9, 35e9, 0.0, vds);
        simulator.setGate(5e9, 30e9, 0, 20e9, 10e9, vgs);
        
        // 执行模拟
        simulator.simulate();
        
        // 获取并保存结果
        double ids = simulator.getDrainCurrent();
        output << vgs << "," << vds << "," << ids << std::endl;
    }
}
```

 4.3 自定义量子参数
可以调整量子修正因子参数来适应不同的物理模型：
```cpp
// 更改量子修正因子参数
simulator.setQuantumFactorParameters(
    4.0,    // 基准修正系数γ₀
    8e9     // 电场依赖系数α
);
```

 5. 结果分析
 5.1 生成的文件格式
模拟器生成的主要结果文件为CSV格式，包含以下列：

 `x(nm),y(nm),z(nm)`  空间坐标
 `potential(V)`  电势
 `electron_density(cm3)`  电子密度
 `hole_density(cm3)`  空穴密度
 `quantum_factor`  量子修正因子
 `electric_field(V/cm)`  电场强度

 5.2 IV特性分析
通过绘制IV特性曲线，可以分析：
1. **阈值电压**：从转移特性曲线确定
2. **亚阈值摆幅**：亚阈值区域的斜率
3. **漏电导**：输出特性曲线中Id对Vds的微分
4. **跨导**：转移特性曲线中Id对Vgs的微分

 5.3 量子效应的影响
通过比较不同量子修正参数的结果，可以分析量子效应对器件性能的影响：
1. 量子限制导致的阈值电压漂移
2. 有效通道厚度的减小
3. 界面附近载流子密度的变化
4. 栅极电容的变化

 6. 高级应用
 6.1 添加新的物理模型
框架设计为模块化结构，可以通过扩展以下方法添加新的物理模型：
 `updateQuantumCorrectionFactor()`  修改量子效应模型
 `updateCarrierDensity()`  使用更复杂的载流子统计
 `updateInterfaceStates()`  添加界面态动力学
 `calculateCurrent()`  实现更精确的电流计算方法

 6.2 网格优化
对于复杂结构，可以实现自适应网格生成：
```cpp
// 生成根据电场梯度优化的非均匀网格
void generateAdaptiveMesh() {
    // 根据电场梯度确定局部网格密度
    // 在关键区域（如界面附近）使用细网格
    // 在变化缓慢区域使用粗网格
}
```

 6.3 并行计算优化
框架已使用OpenMP实现基本并行，可以进一步优化：
```cpp
// 设置OpenMP并行参数
void setParallelParameters(int num_threads) {
    omp_set_num_threads(num_threads);
}
```

 7. 故障排除
 7.1 常见收敛问题
如果模拟无法收敛，尝试以下方法：
1. 减小阻尼因子（默认0.3）
2. 增加最大迭代次数（使用`setMaxIterations()`）
3. 放宽收敛容差（使用`setConvergenceTolerance()`）
4. 使用更细的网格分辨率
5. 减小电压步长，从更低的电压开始迭代

 7.2 数值问题
对于数值计算问题：
1. 检查掺杂浓度是否合理（过高或过低均可能导致问题）
2. 确保电极位置正确且不重叠
3. 边界条件设置是否合理

 8. 后续开发
本框架可以进一步扩展的方向：
1. 添加温度依赖模型
2. 实现全量子非平衡格林函数（NEGF）方法
3. 集成机器学习驱动的自适应网格生成
4. 添加蒙特卡洛模拟模块
5. 实现更完整的界面态动力学模型

 9. 参考文献
1. S. Datta, "Quantum Transport: Atom to Transistor", Cambridge University Press
2. D. Vasileska, S.M. Goodnick, "Computational Electronics", Morgan & Claypool
3. M. Lundstrom, "Fundamentals of Carrier Transport", Cambridge University Press
4. T. Grasser, et al., "Advanced Transport Models for Sub100 nm MOSFETs", IEEE Proc.


   
 编译与运行说明：
 编译要求：
本模拟框架使用C++编写，需要以下依赖：
1. C++14或更高版本
2. Eigen库（版本3.3或更高）
3. OpenMP（用于并行计算）
4. CMake（版本3.10或更高，用于构建工程）

 文件结构
 `QuantumTransportSimulator.h`  主模拟框架头文件
 `devicesetupdemo.cpp`  演示如何设置和模拟单个器件
 `ivcurvegenerator.cpp`  生成IV特性曲线的程序
 `visualization.py`  用于可视化结果的Python脚本

 编译步骤
 1. 创建CMakeLists.txt
在项目根目录创建`CMakeLists.txt`文件，内容如下：
```cmake
cmake_minimum_required(VERSION 3.10)
project(NanoMOSFET_Simulator CXX)

 C++标准设置
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

 启用OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

 Eigen库路径（根据实际情况修改）
include_directories(/usr/include/eigen3)

 添加可执行文件
add_executable(device_demo devicesetupdemo.cpp)
add_executable(iv_generator ivcurvegenerator.cpp)

 链接库
target_link_libraries(device_demo OpenMP::OpenMP_CXX)
target_link_libraries(iv_generator OpenMP::OpenMP_CXX)
```

 2. 构建项目
```bash
 创建并进入构建目录
mkdir build
cd build

 配置
cmake ..

 编译
make
```

 运行示例
 1. 运行单个器件模拟
```bash
./device_demo
```
这将模拟一个纳米尺度MOSFET器件，并生成以下文件：
 `nano_mosfet_results.csv`  包含器件中各点的电势、电荷密度等数据
 `nano_mosfet_results_IV.csv`  包含单个工作点的IV特性

 2. 生成IV特性曲线
```bash
./iv_generator
```
这将在不同栅极电压和漏极电压条件下模拟器件，并生成`iv_curves.csv`文件。

 3. 可视化结果
需要Python及相关包（numpy, matplotlib, pandas）：
```bash
python visualization.py
```
这将生成各种可视化图形，包括：
 `potential_2d_z5.0.png`  z=5nm处的2D电势分布
 `electron_density_z5.0.png`  z=5nm处的电子密度分布
 `quantum_factor_z5.0.png`  z=5nm处的量子修正因子分布
 `potential_along_channel.png`  沿通道方向的电势分布
 `output_characteristics.png`  输出特性曲线（IdVds）
 `transfer_characteristics.png`  转移特性曲线（IdVgs）
 `log_transfer_characteristics.png`  对数比例转移特性曲线

 自定义模拟
您可以通过修改示例文件中的参数来自定义模拟：
1. 修改器件尺寸、材料参数
2. 调整网格分辨率
3. 更改电极位置和电压
4. 设置不同的掺杂分布
5. 修改量子修正因子参数

具体示例请参考`devicesetupdemo.cpp`中的注释说明。


