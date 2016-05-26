import os, sys
import mdtraj as md


def MDTrajParser(fNameTraj):

    skipMDtraj = False
    if fNameTraj and (fNameTraj.endswith("atom") or fNameTraj.endswith("xyz") or fNameTraj.endswith("lammpstrj")):  ## these trajectory styles are handled by LAMMPSParserTraj.py
        skipMDtraj = True

    if fNameTraj and skipMDtraj == False:

        ### LOADING TRAJECTORY AND TOPOLOGY
        mdTrajectory =  md.load(os.path.dirname(os.path.abspath(sys.argv[1])) + '/' + fNameTraj, top='top.pdb')
        mdTopology = md.load_topology('top.pdb')
        ###

        MDTrajAtomPosition = mdTrajectory.xyz                  # atomic positions per frame
        MDTrajAtomPosition = MDTrajAtomPosition.tolist()

        MDTrajSimulationCell = mdTrajectory.unitcell_vectors   # simulation cell per frame
        MDTrajSimulationCell = MDTrajSimulationCell.tolist()

    else:
        MDTrajAtomPosition   = None
        MDTrajSimulationCell = None

    os.remove('top.pdb')  # remove the temporary .pdb topology (required to load a traj file through MDTraj)

    return (MDTrajAtomPosition, MDTrajSimulationCell)