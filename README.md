# This project is meant as a temporary fix to issue relating to issue #399 until a full fix is implemented in trunk.
The VertexBoundaryRefinement modifier was provided by jmosborne and DGermano8 as a possible fix to this issue. 

An example test is given in the test folder. 

This can be left within the project if one wishes to execute this test through the project after the project is correctly built. Otherwise this will need to be place within "Chaste/cell_based/src/simulation/modifiers" if you intend to work from trunk. However, doing this is not recomended.

Within src/ReplaceTheCorrespondingFileInTrunk is a slightly modified version of the MutableVertexMesh file which will allow for new nodes at the boundary to inherit the diffusion radius of their parent node to prevent some errors from occuring.
