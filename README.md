# Structural-FE-codes
(For info about this project go to https://nf55-55.github.io/projects/Efficient%20structural%20FE%20code/)

##Current workflow:

- create the mesh in either: Abaqus, Ansys or SALOME.
	In Abaqus it's possible to write out the .inp file containing also all the mesh info even if the number of nodes/element exceeds the license limit.
- create the sets/groups of nodes on which the boundary conditions (BC) and loads will be applied (fundamental for large/complex meshes to do this
	step in the mesher itself, as the alternative is either to plot the mesh in MATLAB, which is possible only for small meshes, and take note of
	the node IDs, or smartly select the nodes by using for example find(nodes(:,1)==0) -> finding all nodes at x coordinate 0)
- for Abaqus and Ansys write out the .inp file (which contains both the mesh and sets of nodes). 
	While for SALOME export the mesh as .dat file as well as the mesh groups containing the BC and load nodes as .dat files.
- Add the mesh .dat or .inp file to the .\MESH folder, and rename the file with the name of the software used, e.g.: "Abaqus.inp", 
	"Ansys.inp", "SALOME.dat". The only valid names that should be included in the name of the file are: "Abaqus", "Ansys", "SALOME".
	This allows the program to recognize automatically the type of mesh file.
- For Ansys and Abaqus the node sets should be created manually, by copy-pasting the list of node IDs of the .inp file in new .txt files making sure all IDs
	are on the same line, and added to ".\BC and load node sets" folder. Sometimes these programs recognize patterns in the node IDs
	so it's possible that only the first, last and the increment between consecutive ones is given, in that case use the python script in the sets folder to generate the list of node IDs.
	For SALOME just add the .dat files containing the node groups in the ".\BC and load node sets\SALOME dat node groups\" folder.
- IF Ansys or Abaqus nodes sets were used, then load them as described in the .\main and apply the BC and F as you want.
	IF  SALOME node groups were used, then read the script and the one function used to convert them in the usual .txt format contained in
	".\BC and load node sets\SALOME dat node groups\" and change it accordingly. Then load them as described in the .\main and apply the
	BC and F as you want.
- From now on only modification to the .\main parameters should be applied.
- CHECK and modify the main parameters at the top of the script.
- SET material properties
- RUN .\main
- a .vtk file of the results will be saved in ".\VTK RESULTS" folder.
