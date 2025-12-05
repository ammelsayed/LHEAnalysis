# LHEAnalysis

The code implements a general-purpose tool for converting Les Houches Event (LHE) files into ROOT files containing structured physics objects. The tool is designed for truth-level Monte Carlo studies where each event may contain initial-state, intermediate-state, and final-state particles with varying multiplicities.

The converter reads events sequentially using pylhe, extracts the full particle record (PDG ID, status flag, four-momentum components, color and mother indices, lifetime, and helicity), and categorizes particles into species (Electrons, Muons, Taus, Quarks, Gluons, Photons, Bosons, etc.) through a configurable PDG map.

For each particle species appearing in any given event, the algorithm creates a dedicated TClonesArray of Particle objects, where Particle is a custom C++ class with a ROOT dictionary. All collected particles for a species in that event are sorted in descending order of transverse momentum (pT), guaranteeing deterministic ordering of leading, subleading, etc. objects. The event is then written into a single TTree called LHEAnalysis, whose branches correspond to the particle types appearing in the LHE file.

The tool allows the user to control which particle statuses (initial, intermediate, final) are stored in the output ROOT file. The PDG map is configurable, allowing users to extend the list of recognized particles or assign custom quantum numbers. The framework additionally computes derived kinematic variables (pT, η, ϕ) and stores a built-in TLorentzVector for each particle.

The resulting ROOT files are fully compatible with downstream analysis tools, including traditional TTree loops, ROOT RDataFrame and uproot, enabling user-friendly parton-level analysis.

**Note for intermediate (status=2) particles analysis:**

By default, in MadGraph5, only partiles which are *oneshell* are written in the LHE file. The definition of being on-shell depend of the "bw_cut" parameter of the run_card.
All particles with invariant mass between M - bw_cut * Width and M + bw_cut * Width are consider as on-shell and therefore written in the lhe file (see https://answers.launchpad.net/mg5amcnlo/+faq/2173). 
You can also write the intermediate particles if you handel the decays to **Madspin**.

## Analysis Example 1: Phenomenology of boosted top decays

Below are some results of boosted top decays at partonic level, which shows some important insights of the boosted top phenomology that could be helpful in reconstruction level analysis. The code can be found in the `pyroot_example.py` file. The results below were generated with the MadGrapg5 script given in the data folder, but we used here 1M events instead. The `ttbar.lhe.gz` file uploaded in here only contains 10k events due to file size limits on Github.

<img width="1600" height="3000" alt="Characterstics_1M" src="https://github.com/user-attachments/assets/814ce5fb-4519-4b7d-9287-d1de6f7e5a73" />

For more updates please feel welcomed to visit my website at [ammelsayed.tech](https://ammelsayed.tech/).
