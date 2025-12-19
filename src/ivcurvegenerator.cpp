#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "QuantumTransportSimulator.h"

// 生成I-V特性曲线
void generateIVCurves() {
    std::cout << "===== 生成I-V特性曲线 =====" << std::endl;
    
    // 设备参数
    double length = 35e-9;       // 设备长度：35nm
    double width = 20e-9;        // 设备宽度：20nm
    double height = 10e-9;       // 设备高度：10nm
    double m_effective = 0.25;   // 硅中的有效质量
    double tox = 1.2e-9;         // 氧化层厚度：1.2nm
    double eps_ox = 3.9;         // 二氧化硅的相对介电常数
    
    // 网格设置 - 减少点数以加快模拟速度
    int nx = 50;                 // x方向网格点数
    int ny = 30;                 // y方向网格点数
    int nz = 20;                 // z方向网格点数
    
    // 电极位置
    double source_x_start = 0;
    double source_x_end = 5e-9;      // 源区长度：5nm
    double drain_x_start = 30e-9;
    double drain_x_end = 35e-9;      // 漏区长度：5nm
    double gate_x_start = 5e-9;
    double gate_x_end = 30e-9;      // 栅长：25nm
    double gate_y_start = 0;
    double gate_y_end = width;
    double gate_z_position = height; // 栅极位于顶部
    double substrate_z_position = 0;
    double substrate_voltage = 0.0;  // 衬底电压：0V
    
    // 扫描参数
    std::vector<double> gate_voltages = {0.4, 0.6, 0.8, 1.0, 1.2};
    std::vector<double> drain_voltages = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    double source_voltage = 0.0;     // 源极电压：0V
    
    // 结果文件
    std::ofstream outputFile("iv_curves.csv");
    outputFile << "Vgs(V),Vds(V),Ids(A)" << std::endl;
    
    // 执行扫描
    for (double vgs : gate_voltages) {
        std::cout << "模拟栅压 Vgs = " << vgs << " V" << std::endl;
        
        for (double vds : drain_voltages) {
            std::cout << "  - 漏压 Vds = " << vds << " V" << std::endl;
            
            // 创建模拟器实例
            QuantumTransportSimulator simulator(
                length, width, height, nx, ny, nz, m_effective, tox, eps_ox
            );
            
            // 设置收敛参数 - 对于多次运行，适当减少迭代次数
            simulator.setMaxIterations(30);
            simulator.setConvergenceTolerance(5e-5);
            simulator.setQuantumFactorParameters(3.0, 6e9);
            
            // 设置电极
            simulator.setSourceDrain(
                source_x_start, source_x_end,
                drain_x_start, drain_x_end,
                source_voltage, source_voltage + vds  // 漏极电压相对于源极
            );
            
            simulator.setGate(
                gate_x_start, gate_x_end,
                gate_y_start, gate_y_end,
                gate_z_position, source_voltage + vgs  // 栅极电压相对于源极
            );
            
            simulator.setSubstrate(
                substrate_z_position, substrate_voltage
            );
            
            // 设置掺杂分布
            // 通道区：弱掺杂P型
            simulator.setUniformDoping(
                5e-9, 30e-9,                 // 通道区x范围
                0, width,                    // 全宽
                0, 5e-9,                     // z方向范围
                0.0, 1e17                    // 掺杂：P型1e17 cm^-3
            );
            
            // 源区：重掺杂N型
            simulator.setUniformDoping(
                0, 5e-9,                     // 源区x范围
                0, width,                    // 全宽
                0, height,                   // 全高
                5e19, 0.0                    // 掺杂：N型5e19 cm^-3
            );
            
            // 漏区：重掺杂N型
            simulator.setUniformDoping(
                30e-9, 35e-9,                // 漏区x范围
                0, width,                    // 全宽
                0, height,                   // 全高
                5e19, 0.0                    // 掺杂：N型5e19 cm^-3
            );
            
            // 执行模拟
            simulator.simulate();
            
            // 获取结果
            double ids = simulator.getDrainCurrent();
            
            // 保存结果
            outputFile << vgs << "," << vds << "," << ids << std::endl;
            
            std::cout << "    Ids = " << ids << " A" << std::endl;
        }
    }
    
    outputFile.close();
    std::cout << "I-V特性曲线已保存到 iv_curves.csv" << std::endl;
}

int main() {
    generateIVCurves();
    return 0;
}