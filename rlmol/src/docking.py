import meeko
import rdkit.Chem
import rdkit.Chem.AllChem
import shutil
import subprocess
import os
from xml.dom.minidom import parse
from openbabel import pybel
from vina import Vina
import vina


obabel_exec = shutil.which('obabel')
conda_prefix = os.environ['CONDA_PREFIX']
mgltools_python = os.path.join(conda_prefix, 'bin/python2.7')
prepare_gpf_exec = os.path.join(
    conda_prefix, 'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py')
prepare_dpf_exec = os.path.join(
    conda_prefix, 'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py')
summarize_results_exec = os.path.join(
    conda_prefix, 'MGLToolsPckgs/AutoDockTools/Utilities24/summarize_results4.py')
autodock_gpu_exec = shutil.which('autodock_gpu_32wi')
autodock_cpu_exec = shutil.which('autodock4')
autogrid_exec = shutil.which('autogrid4')
prepare_ligand_exec = shutil.which('prepare_ligand4.py')
prepare_receptor_exec = shutil.which('prepare_receptor4.py')


def prepare_dpf(ligand_name, receptor_name):
    args = [mgltools_python, prepare_dpf_exec,
            "-l", ligand_name+".pdbqt",
            "-r", receptor_name+".pdbqt"]
    subprocess.run(args, capture_output=False)


def dock_cpu(ligand_name, receptor_name):
    args = [autodock_cpu_exec,  # path to autodock executable
            "-p", ligand_name+"_"+receptor_name+".dpf",  # docking parameter file: .dpf
            "-l", ligand_name+".dlg",  # docking log file: .dlg
            ]
    subprocess.run(args, capture_output=False)


def get_binding_affinity_cpu() -> float:
    args = [mgltools_python, summarize_results_exec,
            "-d", ".",  # directory to read
            "-b"]  # best docking info only
    subprocess.run(args, capture_output=False)

    with open('summary_of_results_1.0', 'r') as f:
        line = f.readlines()[1]
        ba = float(line.split(',')[4])
        f.close()
    return ba


def save_mol(m, name):
    # saves an rdkit mol object to disk
    copy = rdkit.Chem.RWMol(m)
    m_hydrogens = rdkit.Chem.AddHs(copy)
    rdkit.Chem.AllChem.EmbedMolecule(m_hydrogens)
    print(rdkit.Chem.MolToMolBlock(m_hydrogens),
          file=open(name+'.mol', 'w+'))


def save_mol_babel(m, name):
    mol = pybel.readstring("smi", rdkit.Chem.MolToSmiles(m))
    mol.addh()
    mol.make3D()
    mol.write("mol2", name+".mol2")


def save_mol_sdf(m, name):
    m_hydrogens = rdkit.Chem.AddHs(m)
    rdkit.Chem.AllChem.EmbedMolecule(m_hydrogens)
    rdkit.Chem.AllChem.UFFOptimizeMolecule(m_hydrogens)
    writer = rdkit.Chem.SDWriter(name+".sdf")
    writer.write(m_hydrogens)
    writer.close()


def mol_to_pdb(name):
    # convert .mol to pdb file
    args = [obabel_exec,
            name+'.mol',
            '-opdb',
            '-O', name+'.pdb',
            '-h']
    subprocess.run(args, capture_output=True)


def mol_to_mol2(name):
    args = [obabel_exec,
            name+'.mol',
            '-omol2',
            '-O', name+'.mol2']
    subprocess.run(args, capture_output=True)


def mol2_to_pdbqt(name):
    args = [obabel_exec,
            name+'.mol2',
            '-opdbqt',
            '-O', name+'.pdbqt',
            '-xh']
    subprocess.run(args, capture_output=True)


def mol_to_pdbqt(name):
    args = [obabel_exec,
            name+'.mol',
            '-opdbqt',
            '-O', name+'.pdbqt',
            '-xh']
    subprocess.run(args, capture_output=True)


def sdf_to_pdb(dir, name):
    in_path = os.path.join(dir, name)
    out_path = os.path.join(dir, name)
    args = [obabel_exec,
            in_path+'.sdf',
            '-opdb',
            '-O', out_path+'.pdb',
            '-h']
    subprocess.run(args)


def sdf_to_pdbqt(name):
    args = [obabel_exec,
            name+'.sdf',
            '-opdbqt',
            '-O', name+'.pdbqt',
            '-xh']
    subprocess.run(args)


def pdb_to_pdbqt(name):
    args = [obabel_exec,
            name+'.pdb',
            '-opdbqt',
            '-O', name+'.pdbqt',
            '-xh']
    subprocess.run(args, capture_output=True)


def prepare_ligand(name):
    mol_to_pdb(name)
    args = [mgltools_python,
            prepare_ligand_exec,
            '-l',
            name+'.pdb']
    subprocess.run(args, capture_output=True)


def prepare_ligand_obabel(name):
    mol_to_mol2(name)
    mol2_to_pdbqt(name)


def prepare_receptor(name):
    args = [mgltools_python,
            prepare_receptor_exec,
            '-r',
            name+'.pdb']
    subprocess.run(args, capture_output=True)


def read_grid_file(file_name):
    grid = {}
    with open(file_name) as f:
        lines = f.readlines()

        grid['spacing'] = float(lines[1].split()[1])

        grid['npts'] = [int(lines[2].split()[1]),
                        int(lines[2].split()[2]),
                        int(lines[2].split()[3])]

        grid['center'] = [float(lines[3].split()[1]),
                          float(lines[3].split()[2]),
                          float(lines[3].split()[3])]
    return grid


def prepare_grid(receptor_name, ligand_name, grid):
    gpf_options = []
    gpf_options += ['-p', 'spacing={:.3f}'.format(grid['spacing'])]
    gpf_options += ['-p', 'npts={:d},{:d},{:d}'.format(*grid['npts'])]
    gpf_options += ['-p',
                    'gridcenter={:.2f},{:.2f},{:.2f}'.format(*grid['center'])]

    args = [mgltools_python,  # python version for which mgltools was written
            prepare_gpf_exec,  # path to prepare gpf executable
            "-r", receptor_name+".pdbqt",  # protein file: .pdbqt
            "-l", ligand_name+".pdbqt",  # ligand file: .pdbqt,
            "-o", receptor_name+".gpf"]  # where to save grid parameter file

    subprocess.run(args + gpf_options)

    args = [autogrid_exec,  # autogrid executable
            "-p", receptor_name+".gpf",  # input grid parameter file
            "-l", receptor_name+".gpf_log"]  # output grid log file name

    subprocess.run(args)


def dock(receptor_name, ligand_name, nruns):
    args = [autodock_gpu_exec,  # path to autodock executable
            "-ffile", receptor_name+".maps.fld",  # protein file: .maps.fld
            "-lfile", ligand_name+".pdbqt",  # ligand file: .pdbqt
            "-nrun", str(nruns),  # number of genetic algorithm runs
            "-resnam", ligand_name+".log",  # name for docking output log
            "-dlgoutput", str(0),  # don't write dlg
            "-gbest", str(0)]  # don't write best pose
    subprocess.run(args, capture_output=True)


def reward(name, goal):
    with parse(f'{name}.log.xml') as xml_doc:
        root = xml_doc.documentElement
        clusters = root.getElementsByTagName('cluster')
        ba = clusters[0].attributes['lowest_binding_energy'].value
        # print("Binding affinity is: "+ba+" kCal/Mol")
    if float(ba) <= goal:
        return 1
    else:
        return -1


def get_binding_affinity(name):
    with parse(f'{name}.log.xml') as xml_doc:
        root = xml_doc.documentElement
        clusters = root.getElementsByTagName('cluster')
        ba = clusters[0].attributes['lowest_binding_energy'].value
        return float(ba)


def evaluate_state(mol, target_name, grid, nruns, gpu: bool) -> float:
    ligand_name = 'ligand_'+str(40)
    save_mol(mol, ligand_name)
    prepare_ligand_obabel(ligand_name)
    prepare_grid(target_name, ligand_name, grid)
    if gpu:
        dock(target_name, ligand_name, nruns)
        return get_binding_affinity(ligand_name)
    else:
        prepare_dpf(ligand_name, target_name)
        dock_cpu(ligand_name, target_name)
        return get_binding_affinity_cpu()


def evaluate_state_vina(
        vina: vina.Vina,
        state: str,
        center: list,
        box_size: list,
        exhaustiveness: int,
        n_poses: int) -> float:
    mol = rdkit.Chem.MolFromSmiles(state)
    protonated_mol = rdkit.Chem.AddHs(mol)
    rdkit.Chem.AllChem.EmbedMolecule(protonated_mol)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_mol)
    pdbqt_string = meeko_prep.write_pdbqt_string()

    vina.set_ligand_from_string(pdbqt_string)
    vina.compute_vina_maps(center=center, box_size=box_size)
    vina.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    _ = vina.score()
    energy_minimized = vina.optimize()
    return energy_minimized[0]
