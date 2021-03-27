## Making

`pyrosetta_help` module is just my module for repetitive tasks. One day I'll tidy it up and distribute it
([its repo](https://github.com/matteoferla/pyrosetta_scripts)).

Start...
```python
import pyrosetta
from pyrosetta_help.init_ops import make_option_string, configure_logger

logger = configure_logger()
extra_options= make_option_string(no_optH=False,
                                  ex1=None,
                                  ex2=None,
                                  #mute='all',
                                  ignore_unrecognized_res=True, # raise error!
                                  load_PDB_components=False,
                                  ignore_waters=False)
pyrosetta.init(extra_options=extra_options)
```

Import. I actually make a electron density restrained fast-relaxed model of PDB:1UBQ some time back, so I am using that.

:TODO: find and copy paste that

## Prep
Let's add an isopeptide bond between 11 and 34, then a disulfide between 24 and 52.

To make life easier, I will add the LINK record by hand.

     LINK         NZ  LYS A  11                 CD  GLU A  34                  1.18  

```python
pose = pyrosetta.rosetta.core.import_pose.pose_from_file('1ubq.pre-iso.pdb')
```

### check
Before starting, I need a diagnostic functions that gets the CONN3 residue to a given residue.
Works for LINK-record–type bond (e.g. isopeptides and covalent ligands) and SS bonds.

```python
def get_isopetide_residue(pose: pyrosetta.Pose, resi:int):
    """
    Given a pose and pose residue index, return the CONN3 connecting residue index or zero.
    """
    residue = pose.residue(resi)
    if residue.n_current_residue_connections() != 3:
        return 0
    return residue.connected_residue_at_resconn(3)

assert get_isopetide_residue(11) == 34

def get_con_xyz(pose, res):
    residue = pose.residue(res)
    if residue.n_current_residue_connections() != 3:
        return float('nan')
    con_atomno = residue.residue_connection(3).atomno()
    #con_atomname = residue.atom_name(con_atomno)
    return residue.xyz(con_atomno)

def get_conn_distance(pose, resi_fore:int, resi_aft: int) -> float: 
    """
    Get the distance between two isopeptide residues connecting atoms,
    given their pose residue indices
    """
    assert get_isopetide_residue(resi_fore) == resi_aft
    fore_xyz = get_con_xyz(resi_fore)
    aft_xyz = get_con_xyz(resi_aft)
    d_vector =fore_xyz - aft_xyz
    return d_vector.norm()
    
assert get_conn_distance(11, 34) < 2. # fails until corrected.
```

I added cysteines, so I could add add constraints all at once.

```python
MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
MutateResidue(target=24, new_res='CYS').apply(pose)
MutateResidue(target=52, new_res='CYS').apply(pose)
```

This means I could set some constrains

```python
from pyrosetta_help.common_ops import get_AtomID
# -----------------------------------------------
HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
cl = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t(2)


cl[1]  = AtomPairConstraint(get_AtomID(pose, 'A', 52, 'SG'),
             get_AtomID('A', 24, 'SG'),
             HarmonicFunc(x0_in=2.05, sd_in=0.3))

cl[2]  = AtomPairConstraint(get_AtomID(pose, 'A', 34, 'CD'),
             get_AtomID('A', 11, 'NZ'),
             HarmonicFunc(x0_in=1.33, sd_in=0.1))

cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
cs.add_constraints(cl)

setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
setup.constraint_set(cs)
setup.apply(pose)
```

Also along the way I checked the progress

Check
```python
import nglview as nv

ngl_selection = ' or '.join(map(str, [1, 11, 34, 24, 52, 20]))
print(ngl_selection)
view = nv.show_rosetta(pose)
view.add_representation('hyperball', selection=ngl_selection)
view  # display Jupyter notebook
```

## Isopeptide

I'll just use fast relax to fix the distances

```python
def neighborhood_of_pair(resi_1, resi_2):
    resi1_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(resi_1)
    resi2_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(resi_2)
    or_sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(resi1_sele, resi2_sele)
    NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
    neigh_sele = NeighborhoodResidueSelector(or_sele, distance=12, include_focus_in_subset=True)
    return neigh_sele

pdb2pose = pose.pdb_info().pdb2pose
neigh_sele = neighborhood_of_pair(pdb2pose(chain='A', res=11),
                                  pdb2pose(chain='A', res=34))
n = neigh_sele.apply(pose)
movemap.set_bb(allow_bb=n)
movemap.set_chi(allow_chi=n)
# ----
stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
scorefxn = pyrosetta.create_score_function('ref2015_cart')
scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 10)
relax.cartesian(True)
relax.minimize_bond_angles(True)
relax.minimize_bond_lengths(True)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.apply(pose)
# ----
scorefxn = pyrosetta.create_score_function('ref2015')
scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.apply(pose)
```

## Disulfide
The same for the cysteines, bringing them closer

```python
neigh_sele = neighborhood_of_pair(24, 52)
scorefxn = pyrosetta.create_score_function('ref2015')
scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 15)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.apply(pose)
```

There are some functions, that do some disulfide operations:
```python
pyrosetta.rosetta.core.conformation.form_disulfide(pose.conformation(), 24, 52)
pose.conformation().detect_disulfides()
v = pyrosetta.rosetta.utility.vector1_std_pair_unsigned_long_unsigned_long_t()
pyrosetta.rosetta.core.conformation.disulfide_bonds(pose.conformation(), v)
```
The vector v is vector1_std_pair_unsigned_long_unsigned_long_t[(24, 52)]

## Non-canonical AA

Phosphorylation is curious because in PDB it is a residue, e.g. `SEP`, 
while in Rosetta it is preferentially a patch.

Adding phosphorylation to serine 20.

```python
MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
MutateResidue(target=20, new_res='SER:phosphorylated').apply(pose)
```

Now for a really non-standard one. Norleucine. 

In the database there is a norleucine.

```python
os.listdir(os.path.join(os.path.split(pyrosetta.__file__)[0], 
                        'database', 
                        'chemical', 
                        'residue_type_sets', 
                        'fa_standard', 
                        'residue_types',
                        'l-ncaa'))
```
But uses a different code (`NLU`) than the PDB (`NLE`), which is odd.
And it requires rotamer data, which I'd have to download from the Rosetta Commons page.
So I made my own.

```python
from rdkit_to_params import Params

p = Params.from_smiles('CCCCC(N*)C(*)=O', name='NLE')
p.PROPERTIES.append('ALIPHATIC')
p.PROPERTIES.append('HYDROPHOBIC')
p.dump('norleucine.params')
```

Pretty sure I could have done

```python
p.add_residuetype(pose)
```
But I simply reloaded it from file (I had to check all worked).

```python
from pyrosetta_help.common_ops import pose_from_file
pose = pose_from_file('1ubq.iso.ss.pdb', ['norleucine.params'])
```

I introduced it:

```python
MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
MutateResidue(target=1, new_res='NLE').apply(pose)
resi1_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(pose.pdb_info().pdb2pose(chain='A', res=1))
NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
neigh_sele = NeighborhoodResidueSelector(resi1_sele, distance=9, include_focus_in_subset=True)
movemap = pyrosetta.MoveMap()
movemap.set_bb(allow_bb=resi1_sele.apply(pose))
movemap.set_chi(allow_chi=neigh_sele.apply(pose))
scorefxn = pyrosetta.create_score_function('ref2015')
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.apply(pose)
```