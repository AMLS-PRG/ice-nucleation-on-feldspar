from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator
import numpy as np
from ase import Atoms
from pymatgen.io.cif import CifWriter
from ase.io import read, write
from ase.build import surface, make_supercell
from ase.neighborlist import NeighborList
from ase.build import rotate
from pymatgen.analysis.local_env import CrystalNN
import os


def expand_slab(structure, expansion_matrix):
    """扩展晶胞"""
    print("Before expansion:", len(structure))
    output = make_supercell(structure, expansion_matrix)
    output.write("expanded_structures.cif")
    return make_supercell(structure, expansion_matrix)
    #return structure * expansion_matrix


def reorient_surface(structure, miller_index, layers=1, vacuum=0):
    """根据 Miller 指数切割晶面并添加真空层"""
    lattice = Lattice(structure.get_cell())  # 从 ASE 获取晶格矩阵
    print("print lattice ********************")
    print(lattice)
    print("print lattice ********************")
    species = [atom.symbol for atom in structure]  # 获取原子种类
    positions = structure.get_positions()  # 获取原子坐标
    
    # 将 ASE 笛卡尔坐标转换为 pymatgen 的分数坐标
    frac_coords = lattice.get_fractional_coords(positions)  # 将笛卡尔坐标转换为分数坐标
    #print(frac_coords) 
    #pymatgen_structure = Structure(lattice, species, positions)  # 创建 pymatgen 结构对象
    pymatgen_structure = Structure(lattice, species, frac_coords)  # 创建 pymatgen 结构对象
    angles = pymatgen_structure.lattice.angles  # 返回 alpha, beta, gamma 的列表
    print("Cell angles (degrees):", angles)
    latticepy = pymatgen_structure.lattice  # 从 pymatgen 获取晶格矩阵
    print(f"lattice information: {latticepy}")
    print('print convert structure *****')
    print(pymatgen_structure)

    # 在这里将 pymatgen 结构对象输出为 CIF 文件
    writer = CifWriter(pymatgen_structure)
    cif_filename = "converted_structure.cif"  # 可以根据需要调整文件名
    writer.write_file(cif_filename)  # 保存为 CIF 格式
    print(f"Converted structure saved as CIF: {cif_filename}")
    
    print(pymatgen_structure.lattice) 
    # 使用 Pymatgen 的 SlabGenerator 进行切割, ref https://pymatgen.org/pymatgen.core.html#pymatgen.core.surface.SlabGenerator
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, center_slab=True, lll_reduce=False, reorient_lattice=False)
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, lll_reduce=False, center_slab=False, reorient_lattice=False, primitive=False)
    slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, reorient_lattice=True, center_slab=False, primitive=False, max_normal_search=None)
    # print(slabgen)
    slab = slabgen.get_slab()
    #print(slabs)

    slab_cif_filename = "reorient_structure.cif"
    slab_writer = CifWriter(slab)
    slab_writer.write_file(slab_cif_filename)

    # 返回所有切割的表面
    return slab

def cut_surface_with_constraints(structure, miller_index, layers=1, vacuum=20):
    """根据 Miller 指数切割晶面并添加真空层"""
    lattice = Lattice(structure.get_cell())  # 从 ASE 获取晶格矩阵
    print("print lattice ********************")
    print(lattice)
    print("print lattice ********************")
    species = [atom.symbol for atom in structure]  # 获取原子种类
    positions = structure.get_positions()  # 获取原子坐标
    
    # 将 ASE 笛卡尔坐标转换为 pymatgen 的分数坐标
    frac_coords = lattice.get_fractional_coords(positions)  # 将笛卡尔坐标转换为分数坐标
    #print(frac_coords) 
    #pymatgen_structure = Structure(lattice, species, positions)  # 创建 pymatgen 结构对象
    pymatgen_structure = Structure(lattice, species, frac_coords)  # 创建 pymatgen 结构对象
    angles = pymatgen_structure.lattice.angles  # 返回 alpha, beta, gamma 的列表
    print("Cell angles (degrees):", angles)
    latticepy = pymatgen_structure.lattice  # 从 pymatgen 获取晶格矩阵
    print(f"lattice information: {latticepy}")
    print('print convert structure *****')
    print(pymatgen_structure)
    '''
    # 在这里将 pymatgen 结构对象输出为 CIF 文件
    writer = CifWriter(pymatgen_structure)
    cif_filename = "converted_structure.cif"  # 可以根据需要调整文件名
    writer.write_file(cif_filename)  # 保存为 CIF 格式
    print(f"Converted structure saved as CIF: {cif_filename}")
    '''    
    print(pymatgen_structure.lattice) 
    # 使用 Pymatgen 的 SlabGenerator 进行切割, ref https://pymatgen.org/pymatgen.core.html#pymatgen.core.surface.SlabGenerator
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, center_slab=True, lll_reduce=False, reorient_lattice=False)
    #slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, lll_reduce=False, center_slab=False, reorient_lattice=False, primitive=False)
    slabgen = SlabGenerator(pymatgen_structure, miller_index, layers, vacuum, reorient_lattice=False, center_slab=False, primitive=False, max_normal_search=None)
    # print(slabgen)
    slabs = slabgen.get_slabs(tol=0.01,ftol=0.01)
    #print(slabs)
    for i, slab in enumerate(slabs):
        print(f"Processing slab {i+1}/{len(slabs)}")
    
        writer = CifWriter(slab)
        cif_filename = f"slabs_{i+1}.cif"  # 可以根据需要调整文件名
        writer.write_file(cif_filename)  # 保存为 CIF 格式

    # 返回所有切割的表面
    return slabs



def process_structure(cif_file, miller_index, expansion_matrix, layers=1, vacuum=20):
    """主流程：切割晶体面，计算成键数，检查终端合理性并保存结果"""
    structure = read(cif_file)
        
    # 打印晶胞信息
    print("Original cell parameters (a, b, c, alpha, beta, gamma):")
    print(structure.get_cell_lengths_and_angles())
        
    # 打印原子数和种类
    print(f"Number of atoms: {len(structure)}")
    print(f"Atomic species: {set(atom.symbol for atom in structure)}")
        
    # 打印空间群信息（如果存在）
    spacegroup = structure.info.get('spacegroup', 'Not specified')
    print(f"Spacegroup: {spacegroup}")
    
    lattice = Lattice(structure.get_cell())  # 从 ASE 获取晶格矩阵
    print(f"lattice information: {lattice}")

    reorient_surface(structure, miller_index)
    
    structure = read('reorient_structure.cif')
    # 扩展晶胞
    expanded_structure = expand_slab(structure, expansion_matrix)



    #structure = read('expanded_structures.cif')
    # 切割表面，获取所有表面
    slabs = cut_surface_with_constraints(
        expanded_structure, [0, 0, 1], layers=layers, vacuum=vacuum
        #structure, [0, 0, 1], layers=layers, vacuum=vacuum
    )
    '''
    # 遍历所有的表面
    for i, slab in enumerate(slabs):
        print(f"Processing slab {i+1}/{len(slabs)}")
    
        writer = CifWriter(slab)
        cif_filename = f"slabs_{miller_index}_{i+1}.cif"  # 可以根据需要调整文件名
        writer.write_file(cif_filename)  # 保存为 CIF 格式
    '''



# 示例调用
cif_file = "./Microcline-MS-1979.cif"  # 替换为实际的 CIF 文件路径
miller_index = [1, 1, 0]  # 替换为实际的 Miller 指数
expansion_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 3]]  # 替换为实际的扩展矩阵

process_structure(cif_file, miller_index, expansion_matrix)

