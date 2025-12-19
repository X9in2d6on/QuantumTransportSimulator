#ifndef QUANTUM_TRANSPORT_SIMULATOR_H
#define QUANTUM_TRANSPORT_SIMULATOR_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <omp.h>

// 定义物理常量
namespace PhysicalConstants {
    const double hbar = 1.0545718e-34;    // 约化普朗克常数 (J·s)
    const double q = 1.602176634e-19;     // 基本电荷 (C)
    const double m0 = 9.1093837e-31;      // 电子静止质量 (kg)
    const double kb = 1.380649e-23;       // 玻尔兹曼常数 (J/K)
    const double eps0 = 8.8541878128e-12; // 真空介电常数 (F/m)
    const double T = 300.0;               // 温度 (K)
}

// 量子输运模拟器类
class QuantumTransportSimulator {
public:
    // 构造函数
    QuantumTransportSimulator(
        double length, 
        double width, 
        double height, 
        int nx, 
        int ny, 
        int nz,
        double effective_mass,
        double oxide_thickness,
        double dielectric_constant
    ) : 
        L(length), 
        W(width), 
        H(height), 
        Nx(nx), 
        Ny(ny), 
        Nz(nz),
        m_eff(effective_mass * PhysicalConstants::m0),
        t_ox(oxide_thickness),
        eps_r(dielectric_constant)
    {
        // 初始化网格
        dx = L / (Nx - 1);
        dy = W / (Ny - 1);
        dz = H / (Nz - 1);
        
        total_points = Nx * Ny * Nz;
	    bc_type.assign(total_points, BCType::None);

        
        // 初始化电势，载流子密度和能带结构
        potential = Eigen::VectorXd::Zero(total_points);
        electron_density = Eigen::VectorXd::Zero(total_points);
        hole_density = Eigen::VectorXd::Zero(total_points);

        Nd = Eigen::VectorXd::Zero(total_points);
        Na = Eigen::VectorXd::Zero(total_points);

        // 初始化量子修正因子
        quantum_correction_factor = Eigen::VectorXd::Ones(total_points);
        
        // 初始化界面态密度
        interface_state_density = Eigen::VectorXd::Zero(total_points);
        
        // 默认设置
        max_iterations = 100;
        convergence_tolerance = 1e-6;
        quantum_factor_gamma0 = 3.0;
        quantum_factor_alpha = 5e9;
        
        std::cout << "量子输运模拟器初始化完成" << std::endl;
        std::cout << "设备尺寸: " << L*1e9 << " nm x " << W*1e9 << " nm x " << H*1e9 << " nm" << std::endl;
        std::cout << "网格分辨率: " << Nx << " x " << Ny << " x " << Nz << " 点" << std::endl;
    }

    // 设置源极和漏极位置与电压
    void setSourceDrain(
        double source_x_start, double source_x_end,
        double drain_x_start, double drain_x_end,
        double source_voltage, double drain_voltage
    ) {
        src_x_start = std::max(0.0, source_x_start);
        src_x_end = std::min(L, source_x_end);
        drn_x_start = std::max(0.0, drain_x_start);
        drn_x_end = std::min(L, drain_x_end);
        
        Vs = source_voltage;
        Vd = drain_voltage;
        
        // 设置源极和漏极的边界条件
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    int idx = getIndex(i, j, k);
                    
                    // 源极区域
                    if (x >= src_x_start && x <= src_x_end) {
                            bc_type[idx] = BCType::Source;
    				potential(idx) = Vs;
                    }
                    // 漏极区域
                    else if (x >= drn_x_start && x <= drn_x_end) {
			bc_type[idx] = BCType::Drain;
                        potential(idx) = Vd;
                    }
                }
            }
        }
        
        std::cout << "设置源极和漏极完成" << std::endl;
        std::cout << "源极位置: " << src_x_start*1e9 << " - " << src_x_end*1e9 << " nm" << std::endl;
        std::cout << "漏极位置: " << drn_x_start*1e9 << " - " << drn_x_end*1e9 << " nm" << std::endl;
        std::cout << "源极电压: " << Vs << " V" << std::endl;
        std::cout << "漏极电压: " << Vd << " V" << std::endl;
    }
    
    // 设置栅极位置与电压
    void setGate(
        double gate_x_start, double gate_x_end,
        double gate_y_start, double gate_y_end,
        double gate_z_position,
        double gate_voltage
    ) {
        this->gate_x_start = std::max(0.0, gate_x_start);
        this->gate_x_end   = std::min(L, gate_x_end);
        this->gate_y_start = std::max(0.0, gate_y_start);
        this->gate_y_end   = std::min(W, gate_y_end);
        this->gate_z_pos   = std::min(H, gate_z_position);

        
        Vg = gate_voltage;
        
        // 设置栅极的边界条件
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            if (x >= gate_x_start && x <= gate_x_end) {
                for (int j = 0; j < Ny; j++) {
                    double y = j * dy;
                    if (x >= this->gate_x_start && x <= this->gate_x_end){
                        // 找到最接近gate_z_pos的网格点
                        int k_gate = static_cast<int>(gate_z_pos / dz + 0.5);
                        k_gate = std::min(k_gate, Nz - 1);
                        
                        int idx = getIndex(i, j, k_gate);
			            bc_type[idx] = BCType::Gate;
			            potential(idx) = Vg;

                    }
                }
            }
        }
        
        std::cout << "设置栅极完成" << std::endl;
        std::cout << "栅极位置: x=" << gate_x_start*1e9 << " - " << gate_x_end*1e9 << " nm, ";
        std::cout << "y=" << gate_y_start*1e9 << " - " << gate_y_end*1e9 << " nm, ";
        std::cout << "z=" << gate_z_pos*1e9 << " nm" << std::endl;
        std::cout << "栅极电压: " << Vg << " V" << std::endl;
    }
    
    // 设置衬底位置与电压
    void setSubstrate(double substrate_z_position, double substrate_voltage) {
        sub_z_pos = std::max(0.0, substrate_z_position);
        Vsub = substrate_voltage;
        
        // 设置衬底的边界条件
        int k_sub = static_cast<int>(sub_z_pos / dz + 0.5);
        k_sub = std::max(0, k_sub);
        
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                int idx = getIndex(i, j, k_sub);
		bc_type[idx] = BCType::Substrate;
                potential(idx) = Vsub;
            }
        }
        
        std::cout << "设置衬底完成" << std::endl;
        std::cout << "衬底位置: z=" << sub_z_pos*1e9 << " nm" << std::endl;
        std::cout << "衬底电压: " << Vsub << " V" << std::endl;
    }
    
    // 设置掺杂分布
    void setDoping(
        const Eigen::VectorXd& donor_concentration,
        const Eigen::VectorXd& acceptor_concentration
    ) {
        if (donor_concentration.size() != total_points || 
            acceptor_concentration.size() != total_points) {
            throw std::invalid_argument("掺杂浓度向量尺寸必须与网格点数匹配");
        }
        
        Nd = donor_concentration;
        Na = acceptor_concentration;
        
        std::cout << "设置掺杂分布完成" << std::endl;
    }
    
    // 设置给定区域的均匀掺杂
    void setUniformDoping(
        double x_start, double x_end,
        double y_start, double y_end,
        double z_start, double z_end,
        double donor_conc, double acceptor_conc
    ) {
        x_start = std::max(0.0, x_start);
        x_end = std::min(L, x_end);
        y_start = std::max(0.0, y_start);
        y_end = std::min(W, y_end);
        z_start = std::max(0.0, z_start);
        z_end = std::min(H, z_end);
        
        // Nd = Eigen::VectorXd::Zero(total_points);
        // Na = Eigen::VectorXd::Zero(total_points);
        
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            if (x >= x_start && x <= x_end) {
                for (int j = 0; j < Ny; j++) {
                    double y = j * dy;
                    if (y >= y_start && y <= y_end) {
                        for (int k = 0; k < Nz; k++) {
                            double z = k * dz;
                            if (z >= z_start && z <= z_end) {
                                int idx = getIndex(i, j, k);
                                Nd(idx) = donor_conc;
                                Na(idx) = acceptor_conc;
                            }
                        }
                    }
                }
            }
        }
        
        std::cout << "设置均匀掺杂完成" << std::endl;
        std::cout << "掺杂区域: x=" << x_start*1e9 << " - " << x_end*1e9 << " nm, ";
        std::cout << "y=" << y_start*1e9 << " - " << y_end*1e9 << " nm, ";
        std::cout << "z=" << z_start*1e9 << " - " << z_end*1e9 << " nm" << std::endl;
        std::cout << "施主浓度: " << donor_conc << " cm^-3" << std::endl;
        std::cout << "受主浓度: " << acceptor_conc << " cm^-3" << std::endl;
    }
    
    // 模拟主函数
    void simulate() {
        std::cout << "开始量子输运模拟..." << std::endl;
        
        // 初始化电势分布
        initializePotential();
	    applyElectrodeBoundaryConditions(); // ★关键：初始化后立刻恢复 Vs/Vd/Vg/Vsub
        
        // 主循环，进行自洽求解
        bool converged = false;
        int iter = 0;
        
        while (!converged && iter < max_iterations) {
            applyElectrodeBoundaryConditions(); // ★关键：每轮开始先钉住边界
       		Eigen::VectorXd previous_potential = potential;
            
            // 1. 计算量子修正因子
            updateQuantumCorrectionFactor();
            
            // 2. 根据电势计算载流子密度
            updateCarrierDensity();
            
            // 3. 考虑界面态影响
            updateInterfaceStates();
            
            // 4. 求解泊松方程
            solvePoissonEquation();
            
	        applyElectrodeBoundaryConditions(); // ★关键：阻尼更新后再钉一次

            // 5. 计算电流
            calculateCurrent();
            
            // 6. 检查收敛性
            double diff_norm = (potential - previous_potential).norm() / potential.norm();
            converged = diff_norm < convergence_tolerance;
            
            std::cout << "迭代 " << iter + 1 << " - 相对变化: " << diff_norm << std::endl;
            
            iter++;
        }
        
        if (converged) {
            std::cout << "模拟成功收敛，共 " << iter << " 次迭代" << std::endl;
        } else {
            std::cout << "警告：达到最大迭代次数 " << max_iterations << "，模拟未收敛" << std::endl;
        }
        
        // 计算最终结果
        calculateFinalResults();
    }
    
    // 保存结果到文件
    void saveResults(const std::string& filename) {
        std::ofstream outfile(filename);
        
        if (!outfile.is_open()) {
            std::cerr << "无法打开文件: " << filename << std::endl;
            return;
        }
        
        // 输出文件头
        outfile << "x(nm),y(nm),z(nm),potential(V),electron_density(cm-3),hole_density(cm-3),";
        outfile << "quantum_factor,electric_field(V/cm)" << std::endl;
        
        // 输出数据
        for (int i = 0; i < Nx; i++) {
            double x = i * dx * 1e9; // 转为nm
            for (int j = 0; j < Ny; j++) {
                double y = j * dy * 1e9; // 转为nm
                for (int k = 0; k < Nz; k++) {
                    double z = k * dz * 1e9; // 转为nm
                    int idx = getIndex(i, j, k);
                    
                    outfile << x << "," << y << "," << z << ","
                            << potential(idx) << ","
                            << electron_density(idx) << ","
                            << hole_density(idx) << ","
                            << quantum_correction_factor(idx) << ","
                            << electric_field(idx) << std::endl;
                }
            }
        }
        
        outfile.close();
        std::cout << "结果已保存到: " << filename << std::endl;
        
        // 保存I-V特性
        std::string iv_filename = filename.substr(0, filename.find_last_of('.')) + "_IV.csv";
        std::ofstream iv_file(iv_filename);
        
        if (!iv_file.is_open()) {
            std::cerr << "无法打开文件: " << iv_filename << std::endl;
            return;
        }
        
        // 输出I-V特性
        iv_file << "Vds(V),Ids(A),Vgs(V)" << std::endl;
        iv_file << Vd - Vs << "," << drain_current << "," << Vg << std::endl;
        
        iv_file.close();
        std::cout << "I-V特性已保存到: " << iv_filename << std::endl;
    }
    
    // 获取模拟结果
    double getDrainCurrent() const {
        return drain_current;
    }
    
    Eigen::VectorXd getPotential() const {
        return potential;
    }
    
    Eigen::VectorXd getElectronDensity() const {
        return electron_density;
    }
    
    // 更改模拟参数
    void setMaxIterations(int max_iter) {
        max_iterations = max_iter;
    }
    
    void setConvergenceTolerance(double tol) {
        convergence_tolerance = tol;
    }
    
    void setQuantumFactorParameters(double gamma0, double alpha) {
        quantum_factor_gamma0 = gamma0;
        quantum_factor_alpha = alpha;
    }

private:
    // 设备物理参数
    double L;                     // 设备长度 (m)
    double W;                     // 设备宽度 (m)
    double H;                     // 设备高度 (m)
    double m_eff;                 // 有效质量 (kg)
    double t_ox;                  // 氧化层厚度 (m)
    double eps_r;                 // 相对介电常数
    
    // 网格参数
    int Nx, Ny, Nz;              // 各方向网格点数
    double dx, dy, dz;           // 网格间距 (m)
    int total_points;            // 总网格点数
    
    // 电极参数
    double src_x_start, src_x_end;  // 源极位置
    double drn_x_start, drn_x_end;  // 漏极位置
    double gate_x_start, gate_x_end; // 栅极X位置
    double gate_y_start, gate_y_end; // 栅极Y位置
    double gate_z_pos;             // 栅极Z位置
    double sub_z_pos;              // 衬底位置
    
    // 电压参数
    double Vs;                    // 源极电压 (V)
    double Vd;                    // 漏极电压 (V)
    double Vg;                    // 栅极电压 (V)
    double Vsub;                  // 衬底电压 (V)

enum class BCType : uint8_t { None=0, Source, Drain, Gate, Substrate };
std::vector<BCType> bc_type;

    
    // 物理量分布
    Eigen::VectorXd potential;           // 电势 (V)
    Eigen::VectorXd electron_density;    // 电子密度 (cm^-3)
    Eigen::VectorXd hole_density;        // 空穴密度 (cm^-3)
    Eigen::VectorXd Nd;                  // 施主浓度 (cm^-3)
    Eigen::VectorXd Na;                  // 受主浓度 (cm^-3)
    Eigen::VectorXd quantum_correction_factor; // 量子修正因子
    Eigen::VectorXd interface_state_density;   // 界面态密度 (cm^-2)
    Eigen::VectorXd electric_field;      // 电场强度 (V/cm)
    
    // 结果
    double drain_current;         // 漏电流 (A)
    
    // 求解参数
    int max_iterations;           // 最大迭代次数
    double convergence_tolerance; // 收敛容差
    double quantum_factor_gamma0; // 量子修正因子参数γ₀
    double quantum_factor_alpha;  // 量子修正因子参数α
    
    // 一维索引 (i,j,k) -> idx 转换
    int getIndex(int i, int j, int k) const {
        return i + j * Nx + k * Nx * Ny;
    }
    
    // 反向索引 idx -> (i,j,k) 转换
    void getCoordinates(int idx, int& i, int& j, int& k) const {
        i = idx % Nx;
        j = (idx / Nx) % Ny;
        k = idx / (Nx * Ny);
    }

    
    // 初始化电势分布
    void initializePotential() {
        // 初始化为线性插值
        for (int i = 0; i < Nx; i++) {
            double x = i * dx;
            double x_normalized = x / L;
            
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    int idx = getIndex(i, j, k);
                    
                    // 线性插值源极到漏极的电势
                    potential(idx) = Vs + (Vd - Vs) * x_normalized;
                }
            }
        }
        

        // 初始化电场强度
        electric_field = Eigen::VectorXd::Zero(total_points);
        
        std::cout << "电势分布初始化完成" << std::endl;
    }
    

    //新增：强制施加电极边界电势
    void applyElectrodeBoundaryConditions() {
        for (int idx = 0; idx < total_points; ++idx) {
            switch (bc_type[idx]) {
                case BCType::Source:    potential(idx) = Vs;   break;
                case BCType::Drain:     potential(idx) = Vd;   break;
                case BCType::Gate:      potential(idx) = Vg;   break;
                case BCType::Substrate: potential(idx) = Vsub; break;
                default: break;
            }
        }
    }

    // 更新量子修正因子
    void updateQuantumCorrectionFactor() {
        // 计算电场强度
        calculateElectricField();
        
        // 根据电场强度更新量子修正因子
        // 公式: γ = γ₀ * exp(α * |∇φ|)
        for (int idx = 0; idx < total_points; idx++) {
            double field_magnitude = electric_field(idx); // V/cm
            field_magnitude *= 100.0; // -> V/m

            double arg = quantum_factor_alpha * field_magnitude;

            // 防止exp溢出：50左右 exp(50) ~ 3e21，已经足够“很大”
            if (arg > 50.0) arg = 50.0;
            if (arg < -50.0) arg = -50.0;

            quantum_correction_factor(idx) = quantum_factor_gamma0 * std::exp(arg);

            // 再保险：避免gamma极端值导致数值问题
            if (quantum_correction_factor(idx) < 1.0)  quantum_correction_factor(idx) = 1.0;
            if (quantum_correction_factor(idx) > 1e6)  quantum_correction_factor(idx) = 1e6;

        }
    }
    
    // 计算电场强度
    void calculateElectricField() {
        // 使用中心差分计算电场
        electric_field = Eigen::VectorXd::Zero(total_points);
        
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    int idx = getIndex(i, j, k);
                    
                    // x方向电场
                    double Ex = (potential(getIndex(i+1, j, k)) - potential(getIndex(i-1, j, k))) / (2 * dx);
                    
                    // y方向电场
                    double Ey = (potential(getIndex(i, j+1, k)) - potential(getIndex(i, j-1, k))) / (2 * dy);
                    
                    // z方向电场
                    double Ez = (potential(getIndex(i, j, k+1)) - potential(getIndex(i, j, k-1))) / (2 * dz);
                    
                    // 电场强度 (单位: V/m)
                    double E_magnitude = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
                    
                    // 转换为V/cm
                    electric_field(idx) = E_magnitude * 0.01;
                }
            }
        }
    }
    
    // 更新载流子密度
    void updateCarrierDensity() {
        // 定义常数
        double kT = PhysicalConstants::kb * PhysicalConstants::T;
        double Nc = 2.0 * std::pow(m_eff * kT / (2 * M_PI * PhysicalConstants::hbar * PhysicalConstants::hbar), 1.5);
        double Nv = Nc; // 简化假设
        double Eg = 1.12 * PhysicalConstants::q; // Si带隙 (J)
        double ni = std::sqrt(Nc * Nv) * std::exp(-Eg / (2 * kT)); // 本征载流子浓度
        
        // 考虑量子效应的载流子密度计算
        #pragma omp parallel for
        for (int idx = 0; idx < total_points; idx++) {
            // 经典玻尔兹曼统计近似
            double phi = potential(idx);
            
            // 考虑量子修正
            double gamma = quantum_correction_factor(idx);
            
            // 计算电子密度 (适用于n型区域)
            electron_density(idx) = ni * std::exp(PhysicalConstants::q * phi / (kT * gamma));
            
            // 计算空穴密度 (适用于p型区域)
            hole_density(idx) = ni * std::exp(-PhysicalConstants::q * phi / (kT * gamma));
            
            // 考虑掺杂的影响
		    int i, j, k;
		    getCoordinates(idx, i, j, k);
		    double x = i * dx;

		    // 只在源/漏接触区使用“掺杂钉死”的近似，避免channel失去栅控
		    bool in_contact = (x >= src_x_start && x <= src_x_end) ||
                     		 (x >= drn_x_start && x <= drn_x_end);

		    if (in_contact && (Nd(idx) > 0 || Na(idx) > 0)) {
   		    	if (Nd(idx) > Na(idx) && Nd(idx) > 0) {
    	    		    electron_density(idx) = (Nd(idx) - Na(idx)) + (ni * ni) / Nd(idx);
       	    			hole_density(idx)     = (ni * ni) / Nd(idx);
   		    	} else if (Na(idx) > 0) {
       	    		    hole_density(idx)     = (Na(idx) - Nd(idx)) + (ni * ni) / Na(idx);
            				electron_density(idx) = (ni * ni) / Na(idx);
   		    	}
		    }//保持前面用phi/gamma算出来的 electron_density / hole_density
        }
    }
    
    // 更新界面态
    void updateInterfaceStates() {
        // 在此实现，考虑Shockley-Read-Hall理论
        // 简化模拟中，对界面态的处理简化为立即达到平衡
        // 实际模型需要考虑界面态的捕获和释放动力学
    }
    
    // 求解泊松方程
    void solvePoissonEquation() {
        // 构建系数矩阵和右侧向量
        Eigen::SparseMatrix<double> A(total_points, total_points);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(total_points);
        
        // 填充系数矩阵
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(7 * total_points); // 每个点最多有7个非零元素
        
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    int idx = getIndex(i, j, k);
                    
                    // 检查边界条件
                    if ((i == 0 || i == Nx-1) || 
                        (j == 0 || j == Ny-1) || 
                        (k == 0 || k == Nz-1) ||
                        // 源极、漏极、栅极和衬底区域
                      (bc_type[idx] != BCType::None))
                    {
                        
                        // 边界点：保持电势固定
                        tripletList.push_back(Eigen::Triplet<double>(idx, idx, 1.0));
                        b(idx) = potential(idx);
                    } else {
                        // 内部点：应用有限差分法离散化拉普拉斯算子
                        double h2x = dx * dx;
                        double h2y = dy * dy;
                        double h2z = dz * dz;
                        
                        // 中心点
                        double coef_center = -2.0/h2x - 2.0/h2y - 2.0/h2z;
                        tripletList.push_back(Eigen::Triplet<double>(idx, idx, coef_center));
                        
                        // x方向相邻点
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i-1, j, k), 1.0/h2x));
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i+1, j, k), 1.0/h2x));
                        
                        // y方向相邻点
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i, j-1, k), 1.0/h2y));
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i, j+1, k), 1.0/h2y));
                        
                        // z方向相邻点
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i, j, k-1), 1.0/h2z));
                        tripletList.push_back(Eigen::Triplet<double>(idx, getIndex(i, j, k+1), 1.0/h2z));
                        
                        // 右侧向量：电荷密度
                        double rho = PhysicalConstants::q * (
                            Nd(idx) - Na(idx) + 
                            hole_density(idx) - electron_density(idx) + 
                            interface_state_density(idx)
                        );
                        
                        // 泊松方程右侧项: ∇²φ = -ρ/ε
                        b(idx) = -rho / (PhysicalConstants::eps0 * eps_r);
                    }
                }
            }
        }
        
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        
        // 求解线性方程组 Ax = b
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd solution = solver.solve(b);
        
        // 更新电势
        double damping_factor = 0.3; // 阻尼因子用于提高稳定性
        potential = potential + damping_factor * (solution - potential);
    }
    
    // 计算电流
    void calculateCurrent() {
        // 计算漏电流，通过整合漏区边界上的电流密度
        double current = 0.0;
        
        // 找到漏区边界
        int i_drain = static_cast<int>(drn_x_start / dx);
        
        // 计算通过该边界的电流
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // 计算电子电流和空穴电流
                int idx = getIndex(i_drain, j, k);
                int idx_prev = getIndex(i_drain-1, j, k);
                
                // 电场在x方向
                double E_x = (potential(idx) - potential(idx_prev)) / dx;
                
                // 考虑漂移和扩散电流
                double n_avg = (electron_density(idx) + electron_density(idx_prev)) / 2.0;
                double p_avg = (hole_density(idx) + hole_density(idx_prev)) / 2.0;
                
                // 电子和空穴的迁移率 (简化模型)
                double mu_n = 1400.0 * 1e-4; // m²/(V·s)
                double mu_p = 450.0 * 1e-4;  // m²/(V·s)
                
                // 漂移电流
                double J_drift_n = PhysicalConstants::q * n_avg * mu_n * E_x;
                double J_drift_p = PhysicalConstants::q * p_avg * mu_p * E_x;
                
                // 扩散电流
                double D_n = mu_n * PhysicalConstants::kb * PhysicalConstants::T / PhysicalConstants::q;
                double D_p = mu_p * PhysicalConstants::kb * PhysicalConstants::T / PhysicalConstants::q;
                
                double dn_dx = (electron_density(idx) - electron_density(idx_prev)) / dx;
                double dp_dx = (hole_density(idx) - hole_density(idx_prev)) / dx;
                
                double J_diff_n = PhysicalConstants::q * D_n * dn_dx;
                double J_diff_p = -PhysicalConstants::q * D_p * dp_dx;
                
                // 总电流密度
                double J_total = J_drift_n + J_drift_p + J_diff_n + J_diff_p;
                
                // 积分计算总电流
                current += J_total * dy * dz;
            }
        }
        
        // 保存漏电流
        drain_current = current;
        std::cout << "计算的漏极电流: " << drain_current << " A" << std::endl;
    }
    
    // 计算最终结果
    void calculateFinalResults() {
        // 在此实现任何需要的后处理
        std::cout << "最终漏极电流: " << drain_current << " A" << std::endl;
        std::cout << "漏极-源极电压: " << Vd - Vs << " V" << std::endl;
        std::cout << "栅极-源极电压: " << Vg - Vs << " V" << std::endl;
    }
};

#endif // QUANTUM_TRANSPORT_SIMULATOR_H