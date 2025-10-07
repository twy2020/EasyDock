# EasyDock - 分子对接自动化工具
[README-EN](https://gitlab.igem.org/2025/software-tools/yau-china/-/blob/main/EasyDock/README.md)


[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Molecular Docking](https://img.shields.io/badge/domain-Molecular%20Docking-orange.svg)](https://en.wikipedia.org/wiki/Molecular_docking)

## 📖 项目简介

EasyDock 是一个本地化的分子对接自动化工具，修改自Free_Cloud_Docking(AutoDock colab)项目，移植并优化了分子对接框设置，将复杂的分子对接流程简化为一条命令。基于 Smina 开发，支持全蛋白覆盖对接和多构象搜索，为科研人员提供便捷、高效的分子对接解决方案。

**原项目作者**：https://github.com/quantaosun   
**作者**: Teng  
**邮箱**: tenwonyun@gmail.com  
**GitHub**: [https://github.com/twy2020](https://github.com/twy2020?tab=repositories)

## ✨ 核心特性

### 🎯 全蛋白覆盖对接
- **自动对接盒计算**: 基于蛋白结构自动确定最优对接区域
- **无需手动选择**: 避免传统对接中活性位点选择的偏差
- **全面探索**: 在整个蛋白表面寻找潜在的结合位点

### 🔬 多构象搜索
- **多重配体构象**: 支持配体多构象同时对接
- **构象生成**: 自动从 SMILES 生成 3D 构象
- **构象分析**: 提供构象能量分布和结构多样性分析

### 🛠️ 自动化工作流
- **一键式操作**: 从 PDB ID 和 SMILES 到完整结果
- **智能预处理**: 自动修复蛋白结构，准备配体
- **格式转换**: 自动处理 PDB、PDBQT、SDF 等格式

### 📊 丰富可视化
- **2D 相互作用图**: 详细的蛋白-配体相互作用分析
- **3D 交互视图**: 基于 Web 的 3D 分子查看器
- **PyMOL 会话**: 专业的结构生物学分析会话文件

![dock](https://gitlab.igem.org/2025/software-tools/yau-china/-/raw/main/EasyDock/pic/dock.png)
![result](https://gitlab.igem.org/2025/software-tools/yau-china/-/raw/main/EasyDock/pic/result.png)

## 🚀 快速开始

### 环境要求
- Linux/macOS/Windows (WSL2 推荐用于 Windows)
- Python 3.9+
- 4GB+ 内存
- 10GB+ 磁盘空间

### 一键安装

```bash
# 克隆项目
git clone https://github.com/twy2020/EasyDock.git
cd EasyDock

# 运行自动安装脚本
chmod +x setup_environment.sh
./setup_environment.sh
```

### 手动安装

```bash
# 创建并激活 conda 环境
conda create -n easydock python=3.9 -y
conda activate easydock

# 安装依赖
conda install -c conda-forge openbabel rdkit mdanalysis openmm pdbfixer -y
pip install -r requirements.txt

# 下载 smina
wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina
chmod +x smina
```

## 📖 使用方法

### 基本用法

```bash
# 激活环境
conda activate easydock

# 运行对接（使用默认参数）
python src/main.py --pdb_id 8WRF --smiles "Cn1c2cccc(=O)c-2nc2ccccc21"
```

### 高级参数

```bash
# 自定义对接参数
python src/main.py \
  --pdb_id 8WRF \
  --smiles "Cn1c2cccc(=O)c-2nc2ccccc21" \
  --ligand_name LIG \
  --work_dir my_docking_results \
  --exhaustiveness 64 \
  --num_modes 200
```

### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--pdb_id` | PDB 数据库标识符 (必需) | - |
| `--smiles` | 配体 SMILES 字符串 (必需) | - |
| `--ligand_name` | PDB 中配体的残基名称 | LIG |
| `--work_dir` | 工作目录 | docking_results |
| `--exhaustiveness` | 搜索强度 (8-128) | 32 |
| `--num_modes` | 生成构象数量 | 100 |
| `--config` | 配置文件路径 | - |

### 配置文件

创建 `config.yaml` 文件进行批量配置：

```yaml
# 输入参数
pdb_id: "8WRF"
smiles: "Cn1c2cccc(=O)c-2nc2ccccc21"
ligand_name: "LIG"

# 工作目录
work_dir: "docking_results"

# 对接参数
exhaustiveness: 64
num_modes: 200
energy_range: 5

# 配体准备
num_conformations: 10
```

使用配置文件运行：
```bash
python src/main.py --config config.yaml
```

## 🔧 工作流程

EasyDock 自动执行以下完整流程：

### 1. 📥 数据获取
- 从 RCSB PDB 下载蛋白结构
- 验证 PDB ID 和文件完整性

### 2. 🧬 蛋白准备
- 分离蛋白和配体
- 结构修复和优化
- 添加氢原子和电荷
- 转换为 PDBQT 格式

### 3. ⚗️ 配体准备
- SMILES 到 3D 结构转换
- 多构象生成和优化
- 格式转换和验证

### 4. 🎯 分子对接
- 自动对接盒计算
- 多构象并行对接
- 能量评分和排序

### 5. 📊 结果分析
- 构象能量分析
- 相互作用分析
- 多种格式输出

### 6. 👁️ 可视化
- 2D 相互作用图
- 3D 交互式视图
- PyMOL 会话文件

## 📁 输出文件

运行完成后，工作目录包含类似：

```
docking_results/
├── 8WRF-receptor.pdb              # 准备的蛋白结构
├── receptor.pdbqt                 # 对接用蛋白
├── small.sdf                      # 配体3D结构
├── small_conformation.sdf         # 多构象配体
├── Dockted.pdb                    # 对接结果
├── Dockted.sdf                    # 转换后的结果
├── complex_prepared.pdb           # 复合物结构
├── 2d_interactions.html           # 2D相互作用图
├── 3d_view.html                   # 3D交互视图
├── 3d_session.pse                 # PyMOL会话文件
└── Docked.log                     # 详细对接日志
```

## 🔍 结果解读

### 对接评分
- **结合能**: 负值表示有利结合 (单位: kcal/mol)
- **RMSD**: 构象间结构差异
- **构象排名**: 按结合能排序

### 相互作用分析
- 氢键、疏水相互作用、π-π 堆积等
- 关键残基识别
- 结合模式分析

## 🐛 故障排除

### 常见问题

**Q: 找不到 smina 可执行文件**  
A: 确保 smina 已下载并位于 PATH 中，或使用完整路径

**Q: RDKit 导入错误**  
A: 使用 conda 安装: `conda install -c conda-forge rdkit`

**Q: PDB 下载失败**  
A: 检查网络连接和 PDB ID 有效性

**Q: 对接过程内存不足**  
A: 减少 `num_modes` 参数或使用更大内存机器

### 日志调试

查看详细日志获取错误信息：
```bash
tail -f docking_results/Docked.log
```

## 📚 技术细节

### 算法核心
- **对接引擎**: Smina (AutoDock Vina 分支)
- **构象生成**: RDKit ETKDG 方法
- **蛋白处理**: PDBFixer + OpenBabel
- **可视化**: Py3Dmol + Prolif + MDAnalysis

### 性能优化
- 自动并行处理
- 内存高效的数据结构
- 增量式结果保存

## 🤝 贡献指南

我们欢迎社区贡献！请遵循以下步骤：

1. Fork 本项目
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

## 📄 许可证

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 🙏 致谢

- **Smina 团队**: 提供优秀的分子对接引擎
- **RDKit 社区**: 强大的化学信息学工具包
- **MDAnalysis**: 分子动力学分析工具
- **PyMOL**: 专业的分子可视化软件

## 📞 支持与联系

如果您遇到问题或有建议：

- 📧 邮箱: tenwonyun@gmail.com
- 🐛 [提交 Issue](https://github.com/twy2020/EasyDock/issues)
- 💬 讨论区: [GitHub Discussions](https://github.com/twy2020/EasyDock/discussions)

## 📊 引用

如果您在研究中使用了 EasyDock，请引用：

```bibtex
@software{easydock2024,
  title = {EasyDock: Automated Molecular Docking Pipeline},
  author = {Teng},
  year = {2024},
  url = {https://github.com/twy2020/EasyDock},
  note = {Local automated molecular docking tool with full protein coverage}
}
```