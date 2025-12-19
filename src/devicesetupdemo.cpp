#include <iostream>
#include <string>
#include "QuantumTransportSimulator.h"

int main() {
    std::cout << "===== 纳米尺度MOSFET器件三维量子输运模拟 =====" << std::endl;
    
    // 设备参数
    double length = 35e-9;       // 设备长度：35nm
    double width = 20e-9;        // 设备宽度：20nm
    double height = 10e-9;       // 设备高度：10nm
    double m_effective = 0.25;   // 硅中的有效质量
    double tox = 1.2e-9;         // 氧化层厚度：1.2nm
    double eps_ox = 3.9;         // 二氧化硅的相对介电常数
    
    // 网格设置
    int nx = 70;                 // x方向网格点数
    int ny = 40;                 // y方向网格点数
    int nz = 30;                 // z方向网格点数
    
    // 创建模拟器实例
    QuantumTransportSimulator simulator(
        length, width, height, nx, ny, nz, m_effective, tox, eps_ox
    );
    
    // 设置收敛参数
    simulator.setMaxIterations(50);
    simulator.setConvergenceTolerance(1e-5);
    simulator.setQuantumFactorParameters(3.0, 6e9);
    
    // 设置源极、漏极
    double source_x_start = 0;
    double source_x_end = 5e-9;      // 源区长度：5nm
    double drain_x_start = 30e-9;
    double drain_x_end = 35e-9;      // 漏区长度：5nm
    double source_voltage = 0.0;     // 源极电压：0V
    double drain_voltage = 0.6;      // 漏极电压：0.6V
    
    simulator.setSourceDrain(
        source_x_start, source_x_end,
        drain_x_start, drain_x_end,
        source_voltage, drain_voltage
    );
    
    // 设置栅极
    double gate_x_start = 5e-9;
    double gate_x_end = 30e-9;      // 栅长：25nm
    double gate_y_start = 0;
    double gate_y_end = width;
    double gate_z_position = height; // 栅极位于顶部
    double gate_voltage = 1.2;       // 栅极电压：1.2V
    
    simulator.setGate(
        gate_x_start, gate_x_end,
        gate_y_start, gate_y_end,
        gate_z_position, gate_voltage
    );
    
    // 设置衬底
    double substrate_z_position = 0;
    double substrate_voltage = 0.0;  // 衬底电压：0V
    
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
    
    // 保存结果
    simulator.saveResults("nano_mosfet_results.csv");
    
    // 输出漏电流
    std::cout << "漏电流: " << simulator.getDrainCurrent() << " A" << std::endl;
    
    return 0;
}