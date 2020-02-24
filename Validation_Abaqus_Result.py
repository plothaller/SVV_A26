#Code to evaluate Abaqus results 
#Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

#Import Data from initial aileron
Nodes_Initial_pos = pd.read_csv('Nodes_Initial_pos.txt', delim_whitespace = True, header=None)
Nodes_Initial_pos.columns = ["Node", "X", "Y", "Z"]

Element_based_on_Nodes = pd.read_csv('Element_based_on_Nodes.txt', delim_whitespace = True, header=None)
Element_based_on_Nodes.columns = ["Element", "Node", "Node", "Node", "Node"]

LE_Nodes = pd.read_csv('LE_Nodes.txt', delim_whitespace = True)

TE_Nodes = pd.read_csv('TE_Nodes.txt', delim_whitespace = True)

Spar_Nodes = pd.read_csv('Spar_Nodes.txt', delim_whitespace = True)

Spar_Elements = pd.read_csv('Spar_Elements.txt', delim_whitespace = True)

Skin_Nodes = pd.read_csv('Skin_Nodes.txt', delim_whitespace = True)

Skin_Elements = pd.read_csv('Skin_Elements.txt', delim_whitespace = True)

Hinge_2_Nodes=[7,   8,  27,  29,  30,  31, 162, 163, 164, 165, 166, 553, 554, 555, 556, 557,
               565, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622,
               623, 624, 625, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685,
               686, 687, 688, 689, 690, 721, 723, 724, 725, 726, 727, 728, 729]

#Import Data for case 1
Displacement_case1 = pd.read_csv('Case1_Disp.txt', delim_whitespace = True, header=None)
Displacement_case1.columns =["Node", "Displacement_Total", "X_Disp", "Y_Disp", "Z_Disp"]

#Region 1 = skin
Stresses_1_case1 = pd.read_csv('Case1_Mises.txt', delim_whitespace = True, header=None)
Stresses_1_case1.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]

#Region 2 = Spar
Stresses_2_case1= pd.read_csv('Case1_Mises_reg2.txt', delim_whitespace = True, header=None)
Stresses_2_case1.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]


#Import Data for case 2
Displacement_case2 = pd.read_csv('Bend_Case_Out_U.txt', delim_whitespace = True, header=None)
Displacement_case2.columns =["Node", "Displacement_Total", "X_Disp", "Y_Disp", "Z_Disp"]

#Region 1 = skin
Stresses_1_case2 = pd.read_csv('Bend_Case_Out_Mises.txt', delim_whitespace = True, header=None)
Stresses_1_case2.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]

#Region 2 = Spar
Stresses_2_case2 = pd.read_csv('Bend_Case_Out_Mises_reg2.txt', delim_whitespace = True, header=None)
Stresses_2_case2.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]


#Import Data for case 3
Displacement_case3 = pd.read_csv('Case3_Disp.txt', delim_whitespace = True, header=None)
Displacement_case3.columns =["Node", "Displacement_Total", "X_Disp", "Y_Disp", "Z_Disp"]

#Region 1 = skin
Stresses_1_case3 = pd.read_csv('Case3_Mises.txt', delim_whitespace = True, header=None)
Stresses_1_case3.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]

#Region 2 = Spar
Stresses_2_case3 = pd.read_csv('Case3_Mises_reg2.txt', delim_whitespace = True, header=None)
Stresses_2_case3.columns =["Element", "Int_Pt", "Mises_1", "Mises_2", "S12_1", "S12_2"]


#Calculate average for stresses between Loc1 and Loc2
#Case1 Region1
del Stresses_1_case1['Int_Pt']
avg_1_case1_Mises=(Stresses_1_case1.Mises_1+Stresses_1_case1.Mises_2)/2
del Stresses_1_case1['Mises_1']
del Stresses_1_case1['Mises_2']
Stresses_1_case1['Mises_avg']=avg_1_case1_Mises
avg_1_case1_S12=(Stresses_1_case1.S12_1+Stresses_1_case1.S12_2)/2
del Stresses_1_case1['S12_1']
del Stresses_1_case1['S12_2']
Stresses_1_case1['S12_avg']=avg_1_case1_S12

#Case1 Region2
del Stresses_2_case1['Int_Pt']
avg_2_case1_Mises=(Stresses_2_case1.Mises_1+Stresses_2_case1.Mises_2)/2
del Stresses_2_case1['Mises_1']
del Stresses_2_case1['Mises_2']
Stresses_2_case1['Mises_avg']=avg_2_case1_Mises
avg_2_case1_S12=(Stresses_2_case1.S12_1+Stresses_2_case1.S12_2)/2
del Stresses_2_case1['S12_1']
del Stresses_2_case1['S12_2']
Stresses_2_case1['S12_avg']=avg_2_case1_S12

#Case2 Region1
del Stresses_1_case2['Int_Pt']
avg_1_case2_Mises=(Stresses_1_case2.Mises_1+Stresses_1_case2.Mises_2)/2
del Stresses_1_case2['Mises_1']
del Stresses_1_case2['Mises_2']
Stresses_1_case2['Mises_avg']=avg_1_case2_Mises
avg_1_case2_S12=(Stresses_1_case2.S12_1+Stresses_1_case2.S12_2)/2
del Stresses_1_case2['S12_1']
del Stresses_1_case2['S12_2']
Stresses_1_case2['S12_avg']=avg_1_case2_S12

#Case2 Region2
del Stresses_2_case2['Int_Pt']
avg_2_case2_Mises=(Stresses_2_case2.Mises_1+Stresses_2_case2.Mises_2)/2
del Stresses_2_case2['Mises_1']
del Stresses_2_case2['Mises_2']
Stresses_2_case2['Mises_avg']=avg_2_case2_Mises
avg_2_case2_S12=(Stresses_2_case2.S12_1+Stresses_2_case2.S12_2)/2
del Stresses_2_case2['S12_1']
del Stresses_2_case2['S12_2']
Stresses_2_case2['S12_avg']=avg_2_case2_S12

#Case3 Region1
del Stresses_1_case3['Int_Pt']
avg_1_case3_Mises=(Stresses_1_case3.Mises_1+Stresses_1_case3.Mises_2)/2
del Stresses_1_case3['Mises_1']
del Stresses_1_case3['Mises_2']
Stresses_1_case3['Mises_avg']=avg_1_case3_Mises
avg_1_case3_S12=(Stresses_1_case3.S12_1+Stresses_1_case3.S12_2)/2
del Stresses_1_case3['S12_1']
del Stresses_1_case3['S12_2']
Stresses_1_case3['S12_avg']=avg_1_case3_S12

#Case3 Region2
del Stresses_2_case3['Int_Pt']
avg_2_case3_Mises=(Stresses_2_case3.Mises_1+Stresses_2_case3.Mises_2)/2
del Stresses_2_case3['Mises_1']
del Stresses_2_case3['Mises_2']
Stresses_2_case3['Mises_avg']=avg_2_case3_Mises
avg_2_case3_S12=(Stresses_2_case3.S12_1+Stresses_2_case3.S12_2)/2
del Stresses_2_case3['S12_1']
del Stresses_2_case3['S12_2']
Stresses_2_case3['S12_avg']=avg_2_case3_S12

#Find maximum stresses:
#Case 1
case1_1_max_misses=np.max(Stresses_1_case1.Mises_avg)
case1_1_max_s12=np.max(Stresses_1_case1.S12_avg)
case1_2_max_misses=np.max(Stresses_2_case1.Mises_avg)
case1_2_max_s12=np.max(Stresses_2_case1.S12_avg)
#Case 2
case2_1_max_misses=np.max(Stresses_1_case2.Mises_avg)
case2_1_max_s12=np.max(Stresses_1_case2.S12_avg)
case2_2_max_misses=np.max(Stresses_2_case2.Mises_avg)
case2_2_max_s12=np.max(Stresses_2_case2.S12_avg)
#Case 3
case3_1_max_misses=np.max(Stresses_1_case3.Mises_avg)
case3_1_max_s12=np.max(Stresses_1_case3.S12_avg)
case3_2_max_misses=np.max(Stresses_2_case3.Mises_avg)
case3_2_max_s12=np.max(Stresses_2_case3.S12_avg)
#Case 1 + Case 3
case13_1_max_misses=np.max(Stresses_1_case3.Mises_avg+Stresses_1_case1.Mises_avg)
case13_1_max_s12=np.max(Stresses_1_case3.S12_avg+Stresses_1_case1.S12_avg)
case13_2_max_misses=np.max(Stresses_2_case3.Mises_avg+Stresses_2_case1.Mises_avg)
case13_2_max_s12=np.max(Stresses_2_case3.S12_avg+Stresses_2_case1.S12_avg)

#Find deflection at the hinge line:
Nodes_Spar_position_initial=Spar_Nodes
Nodes_Spar_position_X_initial=[]
Nodes_Spar_position_Y_initial=[]
Nodes_Spar_position_Z_initial=[]
Nodes_Hinge_position_X_initial=[]
Nodes_Hinge_position_Y_initial=[]
Nodes_Hinge_position_Z_initial=[]
Nodes_in_Hinge_line=[]
for x in range(len(Spar_Nodes)):
    if Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']==0:
        Nodes_in_Hinge_line.append(Spar_Nodes.at[x,'Node'])
        Nodes_Hinge_position_X_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X'])
        Nodes_Hinge_position_Y_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y'])
        Nodes_Hinge_position_Z_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z'])
    Nodes_Spar_position_X_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X'])
    Nodes_Spar_position_Y_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y'])
    Nodes_Spar_position_Z_initial.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z'])
Nodes_Spar_position_initial['X']=Nodes_Spar_position_X_initial
Nodes_Spar_position_initial['Y']=Nodes_Spar_position_Y_initial
Nodes_Spar_position_initial['Z']=Nodes_Spar_position_Z_initial
#Case1
Nodes_Spar_position_case1=Spar_Nodes
Nodes_Spar_position_X_case1=[]
Nodes_Spar_position_Y_case1=[]
Nodes_Spar_position_Z_case1=[]
Nodes_Hinge_position_X_case1=[]
Nodes_Hinge_position_Y_case1=[]
Nodes_Hinge_position_Z_case1=[]
for x in range(len(Spar_Nodes)):
    if Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']==0:
        Nodes_Hinge_position_X_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
        Nodes_Hinge_position_Y_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
        Nodes_Hinge_position_Z_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
    Nodes_Spar_position_X_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
    Nodes_Spar_position_Y_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
    Nodes_Spar_position_Z_case1.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case1.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
Nodes_Spar_position_case1['X']=Nodes_Spar_position_X_case1
Nodes_Spar_position_case1['Y']=Nodes_Spar_position_Y_case1
Nodes_Spar_position_case1['Z']=Nodes_Spar_position_Z_case1
#Case2
Nodes_Spar_position_case2=Spar_Nodes
Nodes_Spar_position_X_case2=[]
Nodes_Spar_position_Y_case2=[]
Nodes_Spar_position_Z_case2=[]
Nodes_Hinge_position_X_case2=[]
Nodes_Hinge_position_Y_case2=[]
Nodes_Hinge_position_Z_case2=[]
for x in range(len(Spar_Nodes)):
    if Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']==0:
        Nodes_Hinge_position_X_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
        Nodes_Hinge_position_Y_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
        Nodes_Hinge_position_Z_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
    Nodes_Spar_position_X_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
    Nodes_Spar_position_Y_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
    Nodes_Spar_position_Z_case2.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case2.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
Nodes_Spar_position_case2['X']=Nodes_Spar_position_X_case2
Nodes_Spar_position_case2['Y']=Nodes_Spar_position_Y_case2
Nodes_Spar_position_case2['Z']=Nodes_Spar_position_Z_case2
#Case3
Nodes_Spar_position_case3=Spar_Nodes
Nodes_Spar_position_X_case3=[]
Nodes_Spar_position_Y_case3=[]
Nodes_Spar_position_Z_case3=[]
Nodes_Hinge_position_X_case3=[]
Nodes_Hinge_position_Y_case3=[]
Nodes_Hinge_position_Z_case3=[]
for x in range(len(Spar_Nodes)):
    if Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']==0:
        Nodes_Hinge_position_X_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
        Nodes_Hinge_position_Y_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
        Nodes_Hinge_position_Z_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
    Nodes_Spar_position_X_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'X']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'X_Disp'])
    Nodes_Spar_position_Y_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Y']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'Y_Disp'])
    Nodes_Spar_position_Z_case3.append(Nodes_Initial_pos.at[Spar_Nodes.at[x,'Node']-1,'Z']+Displacement_case3.at[Spar_Nodes.at[x,'Node']-1,'Z_Disp'])
Nodes_Spar_position_case3['X']=Nodes_Spar_position_X_case3
Nodes_Spar_position_case3['Y']=Nodes_Spar_position_Y_case3
Nodes_Spar_position_case3['Z']=Nodes_Spar_position_Z_case3
#Find reference nodes:
Nodes_hingeLine_hinge2=(list(set(Nodes_in_Hinge_line).intersection(set(Hinge_2_Nodes))))
Node_reference_hinge=Nodes_hingeLine_hinge2[-1]

X_reference_hinge=Nodes_Initial_pos.at[Node_reference_hinge-1,'X']
Y_reference_hinge=Nodes_Initial_pos.at[Node_reference_hinge-1,'Y']
Z_reference_hinge=Nodes_Initial_pos.at[Node_reference_hinge-1,'Z']

ind_tip_max=Nodes_Hinge_position_X_initial.index(np.max(Nodes_Hinge_position_X_initial))
ind_tip_min=Nodes_Hinge_position_X_initial.index(np.min(Nodes_Hinge_position_X_initial))

node_tip_max=Nodes_in_Hinge_line[ind_tip_max]
node_tip_min=Nodes_in_Hinge_line[ind_tip_min]

X_initial_hinge_max=Nodes_Initial_pos.at[node_tip_max-1,'X']
Y_initial_hinge_max=Nodes_Initial_pos.at[node_tip_max-1,'Y']
Z_initial_hinge_max=Nodes_Initial_pos.at[node_tip_max-1,'Z']
X_initial_hinge_min=Nodes_Initial_pos.at[node_tip_min-1,'X']
Y_initial_hinge_min=Nodes_Initial_pos.at[node_tip_min-1,'Y']
Z_initial_hinge_min=Nodes_Initial_pos.at[node_tip_min-1,'Z']
#Case1:
X_case1_hinge_max=Displacement_case1.at[node_tip_max-1,'X_Disp']
Y_case1_hinge_max=Displacement_case1.at[node_tip_max-1,'Y_Disp']
Z_case1_hinge_max=Displacement_case1.at[node_tip_max-1,'Z_Disp']
X_case1_hinge_min=Displacement_case1.at[node_tip_min-1,'X_Disp']
Y_case1_hinge_min=Displacement_case1.at[node_tip_min-1,'Y_Disp']
Z_case1_hinge_min=Displacement_case1.at[node_tip_min-1,'Z_Disp']
#Case2:
X_case2_hinge_max=Displacement_case2.at[node_tip_max-1,'X_Disp']
Y_case2_hinge_max=Displacement_case2.at[node_tip_max-1,'Y_Disp']
Z_case2_hinge_max=Displacement_case2.at[node_tip_max-1,'Z_Disp']
X_case2_hinge_min=Displacement_case2.at[node_tip_min-1,'X_Disp']
Y_case2_hinge_min=Displacement_case2.at[node_tip_min-1,'Y_Disp']
Z_case2_hinge_min=Displacement_case2.at[node_tip_min-1,'Z_Disp']
#Case3:
X_case3_hinge_max=Displacement_case3.at[node_tip_max-1,'X_Disp']
Y_case3_hinge_max=Displacement_case3.at[node_tip_max-1,'Y_Disp']
Z_case3_hinge_max=Displacement_case3.at[node_tip_max-1,'Z_Disp']
X_case3_hinge_min=Displacement_case3.at[node_tip_min-1,'X_Disp']
Y_case3_hinge_min=Displacement_case3.at[node_tip_min-1,'Y_Disp']
Z_case3_hinge_min=Displacement_case3.at[node_tip_min-1,'Z_Disp']

#Find bending in Y, Z Direction of the line for 3 cases:
#Case 1:
bending_case1_Y_max=np.degrees(np.arctan(Y_case1_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case1_hinge_max)))))
bending_case1_Y_min=np.degrees(np.arctan(Y_case1_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case1_hinge_min)))))
bending_case1_Z_max=np.degrees(np.arctan(Z_case1_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case1_hinge_max)))))
bending_case1_Z_min=np.degrees(np.arctan(Z_case1_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case1_hinge_min)))))
#Case 2:
bending_case2_Y_max=np.degrees(np.arctan(Y_case2_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case2_hinge_max)))))
bending_case2_Y_min=np.degrees(np.arctan(Y_case2_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case2_hinge_min)))))
bending_case2_Z_max=np.degrees(np.arctan(Z_case2_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case2_hinge_max)))))
bending_case2_Z_min=np.degrees(np.arctan(Z_case2_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case2_hinge_min)))))
#Case 3:
bending_case3_Y_max=np.degrees(np.arctan(Y_case3_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case3_hinge_max)))))
bending_case3_Y_min=np.degrees(np.arctan(Y_case3_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case3_hinge_min)))))
bending_case3_Z_max=np.degrees(np.arctan(Z_case3_hinge_max/(np.abs(X_reference_hinge-X_initial_hinge_max-np.abs(X_case3_hinge_max)))))
bending_case3_Z_min=np.degrees(np.arctan(Z_case3_hinge_min/(np.abs(X_reference_hinge-X_initial_hinge_min-np.abs(X_case3_hinge_min)))))

#Find twist of the aileron
#Initial Data
#Find nodes LE at min X, center X and max X
X_initial_LE_min=0
X_initial_LE_ref=0
X_initial_LE_max=0
Y_initial_LE_min=0
Y_initial_LE_ref=0
Y_initial_LE_max=0
Z_initial_LE_min=0
Z_initial_LE_ref=0
Z_initial_LE_max=0
ind_node_LE_min=0
ind_node_LE_ref=0
ind_node_LE_max=0
for x in range(len(LE_Nodes.Node)):
    if Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']==X_initial_hinge_min:
        X_initial_LE_min=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']
        Y_initial_LE_min=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_LE_min=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Z']
        ind_node_LE_min=LE_Nodes.at[x,'Node']
    if Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']==X_reference_hinge:
        X_initial_LE_ref=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']
        Y_initial_LE_ref=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_LE_ref=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Z']
        ind_node_LE_ref=LE_Nodes.at[x,'Node']
    if Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']==X_initial_hinge_max:
        X_initial_LE_max=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'X']
        Y_initial_LE_max=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_LE_max=Nodes_Initial_pos.at[LE_Nodes.at[x,'Node']-1,'Z']
        ind_node_LE_max=LE_Nodes.at[x,'Node']
              
#Find nodes TE at min X, center X and max X
X_initial_TE_min=0
X_initial_TE_ref=0
X_initial_TE_max=0
Y_initial_TE_min=0
Y_initial_TE_ref=0
Y_initial_TE_max=0
Z_initial_TE_min=0
Z_initial_TE_ref=0
Z_initial_TE_max=0
ind_node_TE_min=0
ind_node_TE_ref=0
ind_node_TE_max=0
for x in range(len(TE_Nodes.Node)):
    if Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']==X_initial_hinge_min:
        X_initial_TE_min=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']
        Y_initial_TE_min=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_TE_min=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Z']
        ind_node_TE_min=TE_Nodes.at[x,'Node']
    if Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']==X_reference_hinge:
        X_initial_TE_ref=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']
        Y_initial_TE_ref=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_TE_ref=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Z']
        ind_node_TE_ref=TE_Nodes.at[x,'Node']
    if Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']==X_initial_hinge_max:
        X_initial_TE_max=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'X']
        Y_initial_TE_max=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Y']
        Z_initial_TE_max=Nodes_Initial_pos.at[TE_Nodes.at[x,'Node']-1,'Z']
        ind_node_TE_max=TE_Nodes.at[x,'Node']

#Case 1
X_case1_LE_min=X_initial_LE_min+Displacement_case1.at[ind_node_LE_min-1,'X_Disp']
X_case1_LE_ref=X_initial_LE_ref+Displacement_case1.at[ind_node_LE_ref-1,'X_Disp']
X_case1_LE_max=X_initial_LE_max+Displacement_case1.at[ind_node_LE_max-1,'X_Disp']
Y_case1_LE_min=Y_initial_LE_min+Displacement_case1.at[ind_node_LE_min-1,'Y_Disp']
Y_case1_LE_ref=Y_initial_LE_ref+Displacement_case1.at[ind_node_LE_ref-1,'Y_Disp']
Y_case1_LE_max=Y_initial_LE_max+Displacement_case1.at[ind_node_LE_max-1,'Y_Disp']
Z_case1_LE_min=Z_initial_LE_min+Displacement_case1.at[ind_node_LE_min-1,'Z_Disp']
Z_case1_LE_ref=Z_initial_LE_ref+Displacement_case1.at[ind_node_LE_ref-1,'Z_Disp']
Z_case1_LE_max=Z_initial_LE_max+Displacement_case1.at[ind_node_LE_max-1,'Z_Disp']
X_case1_TE_min=X_initial_TE_min+Displacement_case1.at[ind_node_TE_min-1,'X_Disp']
X_case1_TE_ref=X_initial_TE_ref+Displacement_case1.at[ind_node_TE_ref-1,'X_Disp']
X_case1_TE_max=X_initial_TE_max+Displacement_case1.at[ind_node_TE_max-1,'X_Disp']
Y_case1_TE_min=Y_initial_TE_min+Displacement_case1.at[ind_node_TE_min-1,'Y_Disp']
Y_case1_TE_ref=Y_initial_TE_ref+Displacement_case1.at[ind_node_TE_ref-1,'Y_Disp']
Y_case1_TE_max=Y_initial_TE_max+Displacement_case1.at[ind_node_TE_max-1,'Y_Disp']
Z_case1_TE_min=Z_initial_TE_min+Displacement_case1.at[ind_node_TE_min-1,'Z_Disp']
Z_case1_TE_ref=Z_initial_TE_ref+Displacement_case1.at[ind_node_TE_ref-1,'Z_Disp']
Z_case1_TE_max=Z_initial_TE_max+Displacement_case1.at[ind_node_TE_max-1,'Z_Disp']

#Case 2
X_case2_LE_min=X_initial_LE_min+Displacement_case2.at[ind_node_LE_min-1,'X_Disp']
X_case2_LE_ref=X_initial_LE_ref+Displacement_case2.at[ind_node_LE_ref-1,'X_Disp']
X_case2_LE_max=X_initial_LE_max+Displacement_case2.at[ind_node_LE_max-1,'X_Disp']
Y_case2_LE_min=Y_initial_LE_min+Displacement_case2.at[ind_node_LE_min-1,'Y_Disp']
Y_case2_LE_ref=Y_initial_LE_ref+Displacement_case2.at[ind_node_LE_ref-1,'Y_Disp']
Y_case2_LE_max=Y_initial_LE_max+Displacement_case2.at[ind_node_LE_max-1,'Y_Disp']
Z_case2_LE_min=Z_initial_LE_min+Displacement_case2.at[ind_node_LE_min-1,'Z_Disp']
Z_case2_LE_ref=Z_initial_LE_ref+Displacement_case2.at[ind_node_LE_ref-1,'Z_Disp']
Z_case2_LE_max=Z_initial_LE_max+Displacement_case2.at[ind_node_LE_max-1,'Z_Disp']
X_case2_TE_min=X_initial_TE_min+Displacement_case2.at[ind_node_TE_min-1,'X_Disp']
X_case2_TE_ref=X_initial_TE_ref+Displacement_case2.at[ind_node_TE_ref-1,'X_Disp']
X_case2_TE_max=X_initial_TE_max+Displacement_case2.at[ind_node_TE_max-1,'X_Disp']
Y_case2_TE_min=Y_initial_TE_min+Displacement_case2.at[ind_node_TE_min-1,'Y_Disp']
Y_case2_TE_ref=Y_initial_TE_ref+Displacement_case2.at[ind_node_TE_ref-1,'Y_Disp']
Y_case2_TE_max=Y_initial_TE_max+Displacement_case2.at[ind_node_TE_max-1,'Y_Disp']
Z_case2_TE_min=Z_initial_TE_min+Displacement_case2.at[ind_node_TE_min-1,'Z_Disp']
Z_case2_TE_ref=Z_initial_TE_ref+Displacement_case2.at[ind_node_TE_ref-1,'Z_Disp']
Z_case2_TE_max=Z_initial_TE_max+Displacement_case2.at[ind_node_TE_max-1,'Z_Disp']

#Case 3
X_case3_LE_min=X_initial_LE_min+Displacement_case3.at[ind_node_LE_min-1,'X_Disp']
X_case3_LE_ref=X_initial_LE_ref+Displacement_case3.at[ind_node_LE_ref-1,'X_Disp']
X_case3_LE_max=X_initial_LE_max+Displacement_case3.at[ind_node_LE_max-1,'X_Disp']
Y_case3_LE_min=Y_initial_LE_min+Displacement_case3.at[ind_node_LE_min-1,'Y_Disp']
Y_case3_LE_ref=Y_initial_LE_ref+Displacement_case3.at[ind_node_LE_ref-1,'Y_Disp']
Y_case3_LE_max=Y_initial_LE_max+Displacement_case3.at[ind_node_LE_max-1,'Y_Disp']
Z_case3_LE_min=Z_initial_LE_min+Displacement_case3.at[ind_node_LE_min-1,'Z_Disp']
Z_case3_LE_ref=Z_initial_LE_ref+Displacement_case3.at[ind_node_LE_ref-1,'Z_Disp']
Z_case3_LE_max=Z_initial_LE_max+Displacement_case3.at[ind_node_LE_max-1,'Z_Disp']
X_case3_TE_min=X_initial_TE_min+Displacement_case3.at[ind_node_TE_min-1,'X_Disp']
X_case3_TE_ref=X_initial_TE_ref+Displacement_case3.at[ind_node_TE_ref-1,'X_Disp']
X_case3_TE_max=X_initial_TE_max+Displacement_case3.at[ind_node_TE_max-1,'X_Disp']
Y_case3_TE_min=Y_initial_TE_min+Displacement_case3.at[ind_node_TE_min-1,'Y_Disp']
Y_case3_TE_ref=Y_initial_TE_ref+Displacement_case3.at[ind_node_TE_ref-1,'Y_Disp']
Y_case3_TE_max=Y_initial_TE_max+Displacement_case3.at[ind_node_TE_max-1,'Y_Disp']
Z_case3_TE_min=Z_initial_TE_min+Displacement_case3.at[ind_node_TE_min-1,'Z_Disp']
Z_case3_TE_ref=Z_initial_TE_ref+Displacement_case3.at[ind_node_TE_ref-1,'Z_Disp']
Z_case3_TE_max=Z_initial_TE_max+Displacement_case3.at[ind_node_TE_max-1,'Z_Disp']

#Calculate twists for 3 cases at both max and min x location wrt reference line
vect_case1_min=np.array([X_case1_TE_min-X_case1_LE_min,Y_case1_TE_min-Y_case1_LE_min,Z_case1_TE_min-Z_case1_LE_min])
vect_case1_ref=np.array([X_case1_TE_ref-X_case1_LE_ref,Y_case1_TE_ref-Y_case1_LE_ref,Z_case1_TE_ref-Z_case1_LE_ref])
vect_case1_max=np.array([X_case1_TE_max-X_case1_LE_max,Y_case1_TE_max-Y_case1_LE_max,Z_case1_TE_max-Z_case1_LE_max])
vect_case2_min=np.array([X_case2_TE_min-X_case2_LE_min,Y_case2_TE_min-Y_case2_LE_min,Z_case2_TE_min-Z_case2_LE_min])
vect_case2_ref=np.array([X_case2_TE_ref-X_case2_LE_ref,Y_case2_TE_ref-Y_case2_LE_ref,Z_case2_TE_ref-Z_case2_LE_ref])
vect_case2_max=np.array([X_case2_TE_max-X_case2_LE_max,Y_case2_TE_max-Y_case2_LE_max,Z_case2_TE_max-Z_case2_LE_max])
vect_case3_min=np.array([X_case3_TE_min-X_case3_LE_min,Y_case3_TE_min-Y_case3_LE_min,Z_case3_TE_min-Z_case3_LE_min])
vect_case3_ref=np.array([X_case3_TE_ref-X_case3_LE_ref,Y_case3_TE_ref-Y_case3_LE_ref,Z_case3_TE_ref-Z_case3_LE_ref])
vect_case3_max=np.array([X_case3_TE_max-X_case3_LE_max,Y_case3_TE_max-Y_case3_LE_max,Z_case3_TE_max-Z_case3_LE_max])
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
#Case 1
twist_case1_min=np.degrees(angle(vect_case1_ref, vect_case1_min))
twist_case1_max=np.degrees(angle(vect_case1_ref, vect_case1_max))
twist_case1_min=np.degrees(np.arctan(vect_case1_min[1]/vect_case1_min[2]))
twist_case1_max=np.degrees(np.arctan(vect_case1_max[1]/vect_case1_max[2]))
#Case 2
twist_case2_min=np.degrees(angle(vect_case2_ref, vect_case2_min))
twist_case2_max=np.degrees(angle(vect_case2_ref, vect_case2_max))
twist_case2_min=np.degrees(np.arctan(vect_case2_min[1]/vect_case2_min[2]))
twist_case2_max=np.degrees(np.arctan(vect_case2_max[1]/vect_case2_max[2]))
#Case 3
twist_case3_min=np.degrees(angle(vect_case3_ref, vect_case3_min))
twist_case3_max=np.degrees(angle(vect_case3_ref, vect_case3_max))
twist_case3_min=np.degrees(np.arctan(vect_case3_min[1]/vect_case3_min[2]))
twist_case3_max=np.degrees(np.arctan(vect_case3_max[1]/vect_case3_max[2]))

#All plots:
#plot original nodes
fig = plt.figure()
initial_plot = plt.axes(projection='3d')
plt.ylim(-300,300 )
initial_plot.set_xlabel('X Label')
initial_plot.set_ylabel('Y Label')
initial_plot.set_zlabel('Z Label')
zdata = np.array(Nodes_Initial_pos.Z)
xdata = np.array(Nodes_Initial_pos.X)
ydata = np.array(Nodes_Initial_pos.Y)
plt.title('Initial Position Nodes')
initial_plot.scatter3D(xdata, ydata, zdata);
#plot displaced nodes case 1
fig = plt.figure()
case1_plt = plt.axes(projection='3d')
plt.ylim(-300,300 )
case1_plt.set_xlabel('X Label')
case1_plt.set_ylabel('Y Label')
case1_plt.set_zlabel('Z Label')
plt.title('Case 1 Position Nodes')
zdata1 = np.array(Displacement_case1.Z_Disp+Nodes_Initial_pos.Z)
xdata1 = np.array(Displacement_case1.X_Disp+Nodes_Initial_pos.X)
ydata1 = np.array(Displacement_case1.Y_Disp+Nodes_Initial_pos.Y)
case1_plt.scatter3D(xdata1, ydata1, zdata1);
#plot displaced nodes case 2
fig = plt.figure()
case2_plt = plt.axes(projection='3d')
plt.ylim(-300,300 )
case2_plt.set_xlabel('X Label')
case2_plt.set_ylabel('Y Label')
case2_plt.set_zlabel('Z Label')
plt.title('Case 2 Position Nodes')
zdata2 = np.array(Displacement_case2.Z_Disp+Nodes_Initial_pos.Z)
xdata2 = np.array(Displacement_case2.X_Disp+Nodes_Initial_pos.X)
ydata2 = np.array(Displacement_case2.Y_Disp+Nodes_Initial_pos.Y)
case2_plt.scatter3D(xdata2, ydata2, zdata2);
#plot displaced nodes case 3
fig = plt.figure()
case3_plt = plt.axes(projection='3d')
plt.ylim(-300,300 )
case3_plt.set_xlabel('X Label')
case3_plt.set_ylabel('Y Label')
case3_plt.set_zlabel('Z Label')
plt.title('Case 3 Position Nodes')
zdata3 = np.array(Displacement_case3.Z_Disp+Nodes_Initial_pos.Z)
xdata3 = np.array(Displacement_case3.X_Disp+Nodes_Initial_pos.X)
ydata3 = np.array(Displacement_case3.Y_Disp+Nodes_Initial_pos.Y)
case3_plt.scatter3D(xdata3, ydata3, zdata3);
#plot initial nodes position Spar
fig = plt.figure()
case_spar_plt = plt.axes(projection='3d')
#case2_spar_plt.set_xlim3d(-3000,3000)
case_spar_plt.set_ylim3d(-300,300)
case_spar_plt.set_zlim3d(-300,300)
case_spar_plt.set_xlabel('X Label')
case_spar_plt.set_ylabel('Y Label')
case_spar_plt.set_zlabel('Z Label')
plt.title('Initial position nodes Hinge Line')
zdata_spar = np.array(Nodes_Hinge_position_Z_initial)
xdata_spar = np.array(Nodes_Hinge_position_X_initial)
ydata_spar = np.array(Nodes_Hinge_position_Y_initial)
case_spar_plt.scatter3D(xdata_spar, ydata_spar, zdata_spar);
#plot displaced nodes case 1 Spar
fig = plt.figure()
case1_spar_plt = plt.axes(projection='3d')
#case2_spar_plt.set_xlim3d(-3000,3000)
case1_spar_plt.set_ylim3d(-300,300)
case1_spar_plt.set_zlim3d(-300,300)
case1_spar_plt.set_xlabel('X Label')
case1_spar_plt.set_ylabel('Y Label')
case1_spar_plt.set_zlabel('Z Label')
plt.title('Case 1 Hinge Line')
zdata1_spar = np.array(Nodes_Hinge_position_Z_case1)
xdata1_spar = np.array(Nodes_Hinge_position_X_case1)
ydata1_spar = np.array(Nodes_Hinge_position_Y_case1)
case1_spar_plt.scatter3D(xdata1_spar, ydata1_spar, zdata1_spar);
#plot displaced nodes case 2 Spar
fig = plt.figure()
case2_spar_plt = plt.axes(projection='3d')
#case2_spar_plt.set_xlim3d(-3000,3000)
case2_spar_plt.set_ylim3d(-300,300)
case2_spar_plt.set_zlim3d(-300,300)
case2_spar_plt.set_xlabel('X Label')
case2_spar_plt.set_ylabel('Y Label')
case2_spar_plt.set_zlabel('Z Label')
plt.title('Case 2 Hinge Line')
zdata2_spar = np.array(Nodes_Hinge_position_Z_case2)
xdata2_spar = np.array(Nodes_Hinge_position_X_case2)
ydata2_spar = np.array(Nodes_Hinge_position_Y_case2)
case2_spar_plt.scatter3D(xdata2_spar, ydata2_spar, zdata2_spar);
#plot displaced nodes case 3 Spar
fig = plt.figure()
case3_spar_plt = plt.axes(projection='3d')
#case2_spar_plt.set_xlim3d(-3000,3000)
case3_spar_plt.set_ylim3d(-300,300)
case3_spar_plt.set_zlim3d(-300,300)
case3_spar_plt.set_xlabel('X Label')
case3_spar_plt.set_ylabel('Y Label')
case3_spar_plt.set_zlabel('Z Label')
plt.title('Case 3 Hinge Line')
zdata3_spar = np.array(Nodes_Hinge_position_Z_case3)
xdata3_spar = np.array(Nodes_Hinge_position_X_case3)
ydata3_spar = np.array(Nodes_Hinge_position_Y_case3)
case3_spar_plt.scatter3D(xdata3_spar, ydata3_spar, zdata3_spar);
#plot case 1 displacement
fig = plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Nodes_Hinge_position_X_case1, Nodes_Hinge_position_Y_case1,'.')
plt.title('Case 1 Deflections')
plt.ylabel('Y Deflection')
plt.subplot(2, 1, 2)
plt.plot(Nodes_Hinge_position_X_case1, Nodes_Hinge_position_Z_case1,'.')
plt.xlabel('X Position')
plt.ylabel('Z Deflection')
plt.show()
#plot case 2 displacement
fig = plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Nodes_Hinge_position_X_case2, Nodes_Hinge_position_Y_case2,'.')
plt.title('Case 2 Deflections')
plt.ylabel('Y Deflection')
plt.subplot(2, 1, 2)
plt.plot(Nodes_Hinge_position_X_case2, Nodes_Hinge_position_Z_case2,'.')
plt.xlabel('X Position')
plt.ylabel('Z Deflection')
plt.show()
#plot case 3 displacement
fig = plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Nodes_Hinge_position_X_case3, Nodes_Hinge_position_Y_case3,'.')
plt.title('Case 3 Deflections')
plt.ylabel('Y Deflection')
plt.subplot(2, 1, 2)
plt.plot(Nodes_Hinge_position_X_case3, Nodes_Hinge_position_Z_case3,'.')
plt.xlabel('X Position')
plt.ylabel('Z Deflection')
plt.show()
#plot all 3 deflections in the same plot
fig = plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Nodes_Hinge_position_X_case1, Nodes_Hinge_position_Y_case1,'o', label='Case 1')
plt.plot(Nodes_Hinge_position_X_case2, Nodes_Hinge_position_Y_case2,'o', label='Case 2')
plt.plot(Nodes_Hinge_position_X_case3, Nodes_Hinge_position_Y_case3,'o', label='Case 3')
plt.plot(Nodes_Hinge_position_X_initial,Nodes_Hinge_position_Y_initial,'o', label='Initial Pos')
plt.title('Deflections from the 3 Cases')
plt.ylabel('Y Deflection')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(Nodes_Hinge_position_X_case1, Nodes_Hinge_position_Z_case1,'o', label='Case 1')
plt.plot(Nodes_Hinge_position_X_case2, Nodes_Hinge_position_Z_case2,'o', label='Case 2')
plt.plot(Nodes_Hinge_position_X_case3, Nodes_Hinge_position_Z_case3,'o', label='Case 3')
plt.plot(Nodes_Hinge_position_X_initial, Nodes_Hinge_position_Z_initial, 'o', label='Initial Pos')
plt.xlabel('X Position')
plt.ylabel('Z Deflection')
plt.legend()
plt.show()
#plot twist lines:
fig = plt.figure()
twist_ref_plot = fig.gca(projection='3d')
z_min = np.array([Z_initial_LE_min,Z_initial_TE_min])
x_min = np.array([X_initial_LE_min,X_initial_TE_min])
y_min = np.array([Y_initial_LE_min,Y_initial_TE_min])
z_ref = np.array([Z_initial_LE_ref,Z_initial_TE_ref])
x_ref = np.array([X_initial_LE_ref,X_initial_TE_ref])
y_ref = np.array([Y_initial_LE_ref,Y_initial_TE_ref])
z_max = np.array([Z_initial_LE_max,Z_initial_TE_max])
x_max = np.array([X_initial_LE_max,X_initial_TE_max])
y_max = np.array([Y_initial_LE_max,Y_initial_TE_max])
twist_ref_plot.set_ylim3d(-300,300)
twist_ref_plot.set_zlim3d(-500,100)
twist_ref_plot.set_xlabel('X Label')
twist_ref_plot.set_ylabel('Y Label')
twist_ref_plot.set_zlabel('Z Label')
plt.title('Reference Twist')
twist_ref_plot.plot(x_min, y_min, z_min)
twist_ref_plot.plot(x_ref, y_ref, z_ref)
twist_ref_plot.plot(x_max, y_max, z_max)
plt.show()
#Case 1 twist
fig = plt.figure()
twist_case1_plot = fig.gca(projection='3d')
z_min_case1 = np.array([Z_case1_LE_min,Z_case1_TE_min])
x_min_case1 = np.array([X_case1_LE_min,X_case1_TE_min])
y_min_case1 = np.array([Y_case1_LE_min,Y_case1_TE_min])
z_ref_case1 = np.array([Z_case1_LE_ref,Z_case1_TE_ref])
x_ref_case1 = np.array([X_case1_LE_ref,X_case1_TE_ref])
y_ref_case1 = np.array([Y_case1_LE_ref,Y_case1_TE_ref])
z_max_case1 = np.array([Z_case1_LE_max,Z_case1_TE_max])
x_max_case1 = np.array([X_case1_LE_max,X_case1_TE_max])
y_max_case1 = np.array([Y_case1_LE_max,Y_case1_TE_max])
twist_case1_plot.set_ylim3d(-300,300)
twist_case1_plot.set_zlim3d(-500,100)
twist_case1_plot.set_xlabel('X Label')
twist_case1_plot.set_ylabel('Y Label')
twist_case1_plot.set_zlabel('Z Label')
plt.title('Case 1 twist')
twist_case1_plot.plot(x_min_case1, y_min_case1, z_min_case1)
twist_case1_plot.plot(x_ref_case1, y_ref_case1, z_ref_case1)
twist_case1_plot.plot(x_max_case1, y_max_case1, z_max_case1)
plt.show()
#Case 2 twist
fig = plt.figure()
twist_case2_plot = fig.gca(projection='3d')
z_min_case2 = np.array([Z_case2_LE_min,Z_case2_TE_min])
x_min_case2 = np.array([X_case2_LE_min,X_case2_TE_min])
y_min_case2 = np.array([Y_case2_LE_min,Y_case2_TE_min])
z_ref_case2 = np.array([Z_case2_LE_ref,Z_case2_TE_ref])
x_ref_case2 = np.array([X_case2_LE_ref,X_case2_TE_ref])
y_ref_case2 = np.array([Y_case2_LE_ref,Y_case2_TE_ref])
z_max_case2 = np.array([Z_case2_LE_max,Z_case2_TE_max])
x_max_case2 = np.array([X_case2_LE_max,X_case2_TE_max])
y_max_case2 = np.array([Y_case2_LE_max,Y_case2_TE_max])
twist_case2_plot.set_ylim3d(-300,300)
twist_case2_plot.set_zlim3d(-500,100)
twist_case2_plot.set_xlabel('X Label')
twist_case2_plot.set_ylabel('Y Label')
twist_case2_plot.set_zlabel('Z Label')
plt.title('Case 2 twist')
twist_case2_plot.plot(x_min_case2, y_min_case2, z_min_case2)
twist_case2_plot.plot(x_ref_case2, y_ref_case2, z_ref_case2)
twist_case2_plot.plot(x_max_case2, y_max_case2, z_max_case2)
plt.show()
#Case 3 twist
fig = plt.figure()
twist_case3_plot = fig.gca(projection='3d')
z_min_case3 = np.array([Z_case3_LE_min,Z_case3_TE_min])
x_min_case3 = np.array([X_case3_LE_min,X_case3_TE_min])
y_min_case3 = np.array([Y_case3_LE_min,Y_case3_TE_min])
z_ref_case3 = np.array([Z_case3_LE_ref,Z_case3_TE_ref])
x_ref_case3 = np.array([X_case3_LE_ref,X_case3_TE_ref])
y_ref_case3 = np.array([Y_case3_LE_ref,Y_case3_TE_ref])
z_max_case3 = np.array([Z_case3_LE_max,Z_case3_TE_max])
x_max_case3 = np.array([X_case3_LE_max,X_case3_TE_max])
y_max_case3 = np.array([Y_case3_LE_max,Y_case3_TE_max])
twist_case3_plot.set_ylim3d(-300,300)
twist_case3_plot.set_zlim3d(-500,100)
twist_case3_plot.set_xlabel('X Label')
twist_case3_plot.set_ylabel('Y Label')
twist_case3_plot.set_zlabel('Z Label')
plt.title('Case 3 twist')
twist_case3_plot.plot(x_min_case3, y_min_case3, z_min_case3)
twist_case3_plot.plot(x_ref_case3, y_ref_case3, z_ref_case3)
twist_case3_plot.plot(x_max_case3, y_max_case3, z_max_case3)
plt.show()

#Outputs
#Print max stresses for all cases
print('Case 1 Maximum Streeses:')
print('Case 1 Region 1 max Mises: ', np.around(case1_1_max_misses,decimals=4), 'Pa')
print('Case 1 Region 2 max Mises: ', np.around(case1_2_max_misses,decimals=4), 'Pa')
print('Case 1 Region 1 max S12: ', np.around(case1_1_max_s12,decimals=4), 'Pa')
print('Case 1 Region 2 max S12: ', np.around(case1_2_max_s12,decimals=4), 'Pa')
print()
print('Case 2 Maximum Streeses:')
print('Case 2 Region 1 max Mises: ', np.around(case2_1_max_misses,decimals=4), 'Pa')
print('Case 2 Region 2 max Mises: ', np.around(case2_2_max_misses,decimals=4), 'Pa')
print('Case 2 Region 1 max S12: ', np.around(case2_1_max_s12,decimals=4), 'Pa')
print('Case 2 Region 2 max S12: ', np.around(case2_2_max_s12,decimals=4), 'Pa')
print()
print('Case 3 Maximum Streeses:')
print('Case 3 Region 1 max Mises: ', np.around(case3_1_max_misses,decimals=4), 'Pa')
print('Case 3 Region 2 max Mises: ', np.around(case3_2_max_misses,decimals=4), 'Pa')
print('Case 3 Region 1 max S12: ', np.around(case3_1_max_s12,decimals=4), 'Pa')
print('Case 3 Region 2 max S12: ', np.around(case3_2_max_s12,decimals=4), 'Pa')
print()
print('Case 1 + Case 3 Maximum Streeses:')
print('Case 1 + Case 3 Region 1 max Mises: ', np.around(case13_1_max_misses,decimals=4), 'Pa')
print('Case 1 + Case 3 Region 2 max Mises: ', np.around(case13_2_max_misses,decimals=4), 'Pa')
print('Case 1 + Case 3 Region 1 max S12: ', np.around(case13_1_max_s12,decimals=4), 'Pa')
print('Case 1 + Case 3 Region 2 max S12: ', np.around(case13_2_max_s12,decimals=4), 'Pa')
print()
#Print hinge line deflection for all cases
print('Bending Results Case 1:')
print('Bending in Y direction for max X: ', np.around(bending_case1_Y_max,decimals=3), 'deg')
print('Bending in Y direction for min X: ', np.around(bending_case1_Y_min,decimals=3), 'deg')
print('Bending in Z direction for max X: ', np.around(bending_case1_Z_max,decimals=3), 'deg')
print('Bending in Z direction for min X: ', np.around(bending_case1_Z_min,decimals=3), 'deg')
print()
print('Bending Results Case 2:')
print('Bending in Y direction for max X: ',np.around(bending_case2_Y_max,decimals=3), 'deg')
print('Bending in Y direction for min X: ',np.around(bending_case2_Y_min,decimals=3), 'deg')
print('Bending in Z direction for max X: ',np.around(bending_case2_Z_max,decimals=3), 'deg')
print('Bending in Z direction for min X: ',np.around(bending_case2_Z_min,decimals=3), 'deg')
print()
print('Bending Results Case 3:')
print('Bending in Y direction for max X: ',np.around(bending_case3_Y_max,decimals=3), 'deg')
print('Bending in Y direction for min X: ',np.around(bending_case3_Y_min,decimals=3), 'deg')
print('Bending in Z direction for max X: ',np.around(bending_case3_Z_max,decimals=3), 'deg')
print('Bending in Z direction for min X: ',np.around(bending_case3_Z_min,decimals=3), 'deg')
print()
print('Bending Results Case 1 + Case 3:')
print('Bending in Y direction for max X: ',np.around(bending_case3_Y_max+bending_case1_Y_max,decimals=3), 'deg')
print('Bending in Y direction for min X: ',np.around(bending_case3_Y_min+bending_case1_Y_min,decimals=3), 'deg')
print('Bending in Z direction for max X: ',np.around(bending_case3_Z_max+bending_case1_Z_max,decimals=3), 'deg')
print('Bending in Z direction for min X: ',np.around(bending_case3_Z_min+bending_case1_Z_min,decimals=3), 'deg')
print()
#Print twist for all cases
print('Twist Results Case 1:')
print('Twist for max X: ', np.around(twist_case1_max,decimals=3), 'deg')
print('Twist for min X: ', np.around(twist_case1_min,decimals=3), 'deg')
print()
print('Twist Results Case 2:')
print('Twist for max X: ',np.around(twist_case2_max,decimals=3), 'deg')
print('Twist for min X: ',np.around(twist_case2_min,decimals=3), 'deg')
print()
print('Twist Results Case 3:')
print('Twist for max X: ',np.around(twist_case3_max,decimals=3), 'deg')
print('Twist for min X: ',np.around(twist_case3_min,decimals=3), 'deg')
print()
print('Twist Results Case 1 + Case 3:')
print('Twist for max X: ',np.around(twist_case1_max+twist_case3_max,decimals=3), 'deg')
print('Twist for min X: ',np.around(twist_case1_min+twist_case3_min,decimals=3), 'deg')
print()