# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""
import numpy as np
import calfem.core as cfc
import triangles as tri
import quads as quad
import calfem.vis as cfv
import coordtransform as ct
import matplotlib.pyplot as plt


def getXY():
    print("Koordinater til X1-X4, med klokken fra nederst i venste hjørne")
    X = np.array(np.zeros(4))
    X[0] = float(input())
    X[1] = float(input())
    X[2] = float(input())
    X[3] = float(input())
    Y = np.array(np.zeros(4))
    Y[0] = float(input("Samme med Y koordinater\n"))
    Y[1] = float(input())
    Y[2] = float(input())
    Y[3] = float(input())

    t = float(input("Tykkelse\n"))
    n = int(input("Elementtype:\n3 Noders trekant: 3\n6 Noders trekant: 6\n4 Noders Firkant: 4\n9 Noders Firkant: 9\n"))

    nx = int(input("Antall noder x-rettning: "))
    ny = int(input("Antall noder y-rettning: "))

    eq = np.array(np.zeros(2))
    f = np.array(np.zeros(2))
    print("Volum last, x så y:")

    eq[0] = float(input())
    eq[1] = float(input())

    print("Last på andre side, x så y:")
    f[0] = float(input())
    f[1] = float(input())

    return X.T, Y.T, t, n, nx, ny, eq, f


#Element Type

# nu = 0.4
# E =5
# Dmat = np.mat([
#         [ 1.0,  nu,  0.],
#         [  nu, 1.0,  0.],
#         [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)
#
# ex = np.array([0,0,2])
# ey = np.array([1,-1,-1])
#
# val = tri.tri6_Kmatrix(ex,ey,Dmat, 0.5)
element = [3,4,6,9]
#X,Y, thickness, numElementNodes, numNodesX, numNodesY,eq, endLoadXY = getXY()
dispx = []
dispy = []
nodesx = []
nodesy = []
elementtype = []

thickness =  0.2

X = np.array([0,8,8,0]).T
Y = np.array([0,0,1,1]).T
eq = np.array([100000,0])
endLoadXY = np.array([0.0,0.0])

#numElementNodes = 9  # Valid numbers 3, 4, 6, 9
for n in element:
    for Nx in range(3, 50, 2):
        numNodesX = Nx
        numNodesY = Nx
        numElementNodes = n

        elTypeInfo= [-1,'Unknown elementtype']

        meshText = 'Unknown mesh'
        if numElementNodes == 3:
            elTypeInfo= [2,'3 Node Triangle']
        elif numElementNodes == 4:
            elTypeInfo= [3,'4 node Quad mesh']
        elif numElementNodes == 6:
            elTypeInfo= [9,'6 node Triangle mesh']
        elif numElementNodes == 9:
            elTypeInfo= [10,'9 node Quad mesh']

        # Number of nodes: Should be odd numbers in order to handle

        #numNodesX = 40
        #numNodesY = 40

        # number of patches that will fit a 9 node element
        numPatchX = (numNodesX-1) // 2
        numPatchX = 1 if numPatchX < 1 else numPatchX
        numPatchY = (numNodesY-1) // 2
        numPatchY = 1 if numPatchY < 1 else numPatchY

        numNodesX = numPatchX*2 + 1
        numNodesY = numPatchY*2 + 1

        if numElementNodes == 6 or numElementNodes == 9:
            numElementsX = (numNodesX-1) // 2
            numElementsY = (numNodesY-1) // 2
        else:
            numElementsX = numNodesX -1
            numElementsY = numNodesY -1

        bDrawMesh = False

        # Cantilever with dimensions H x L x thickness
        #H         =  2.0
        #L         = 10.0






        lx = 1 / (numNodesX-1)
        ly = 1 / (numNodesY-1)

        tmp1 = np.array([[1, X[0], Y[0]],
                             [1, X[1], Y[1]],
                             [1, X[2], Y[2]]])

        tmp2 = np.array([[1, X[2], Y[2]],
                             [1, X[3], Y[3]],
                             [1, X[0], Y[0]]])

        tmp1 = np.linalg.det(tmp1)
        tmp2 = np.linalg.det(tmp2)

        area  = tmp1/2 + tmp2 /2
        # Distributed load in x and y, load pr unit area
        #eq = np.array([0,0])
        #End load, Given as resultant

        #endLoadXY = np.array([1000000.0,0.0])
        #endLoadXY = np.array([3.0e6,0])
        #endLoadXY = np.array([4.2e9,0.0]) # Should give unit disp at Poisson = 0


        eqTotal = eq *area* thickness #Total load for plotting purpose

        # Material properties and thickness

        ep = [1,thickness]
        E  = 2.1e11
        nu = 0.3
        Dmat = np.mat([
                [ 1.0,  nu,  0.],
                [  nu, 1.0,  0.],
                [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

        numNodes    = numNodesX * numNodesY
        numElements = numElementsX * numElementsY
        if numElementNodes in [3,6]:
            numElements *= 2

        #L_elx = L / (numNodesX-1)
        #L_ely = H / (numNodesY-1)

        nelnod = 6

        coords = np.zeros((numNodes,2))
        dofs   = np.zeros((numNodes,2),int)       #Dofs is starting on 1 on first dof
        edofs  = np.zeros((numElements,numElementNodes*2),int) #edofs are also starting on 1 based dof


        inod = 0 # The coords table starts numbering on 0
        idof = 1 # The values of the dofs start on 1
        ndofs = numNodes * 2

        # Set the node coordinates and node dofs




        for i in range(numNodesX):
            for j in range(numNodesY):
                #coords[inod,0] = L_elx * i
                #coords[inod,1] = L_ely * j

                coords[inod, 0] = ct.getCordX(X, Y, lx * i, ly * j)
                coords[inod, 1] = ct.getCordY(X, Y, lx * i, ly * j)


                dofs[inod,0] = idof
                dofs[inod,1] = idof+1
                idof += 2
                inod += 1

        # Set the element connectivites and element dofs
        elnods = np.zeros((numElements,numElementNodes),int)
        eldofs = np.zeros((numElements,numElementNodes*2),int)

        iel = 0
        for ip in range(numPatchX):
            ii = ip*2
            for jp in range(numPatchY):
                jj = jp*2
                # 0 based node numbers, 9 nodes of a 3x3 patch
                nod9 = np.array([
                    (ii  )*numNodesY + (jj  ),
                    (ii+1)*numNodesY + (jj  ),
                    (ii+2)*numNodesY + (jj  ),
                    (ii  )*numNodesY + (jj+1),
                    (ii+1)*numNodesY + (jj+1),
                    (ii+2)*numNodesY + (jj+1),
                    (ii  )*numNodesY + (jj+2),
                    (ii+1)*numNodesY + (jj+2),
                    (ii+2)*numNodesY + (jj+2)],'i')

                if numElementNodes == 3:
                    for i in range(2):
                        for j in range(2):
                            elnods[iel,:] = [nod9[3*i+j],nod9[3*i+j+1],nod9[3*(i+1)+j+1]]
                            iel += 1
                            elnods[iel,:] = [nod9[3*(i+1)+j+1],nod9[3*(i+1)+j],nod9[3*i+j]]
                            iel += 1
                elif numElementNodes == 6:
                    elnods[iel,:] = [nod9[0],nod9[2],nod9[8],nod9[1],nod9[5],nod9[4]]
                    iel += 1
                    elnods[iel,:] = [nod9[8],nod9[6],nod9[0],nod9[7],nod9[3],nod9[4]]
                    iel += 1
                elif numElementNodes == 4:
                    for i in range(2):
                        for j in range(2):
                            elnods[iel,:] = [nod9[3*i+j],nod9[3*i+j+1],nod9[3*(i+1)+j+1],nod9[3*(i+1)+j]]
                            iel += 1
                elif numElementNodes == 9:
                    elnods[iel,:] = [nod9[0],nod9[2],nod9[8],nod9[6],
                                     nod9[1],nod9[5],nod9[7],nod9[3],
                                     nod9[4]]
                    iel += 1


        for iel in range(elnods.shape[0]):
            eldofs[iel, ::2] = elnods[iel,:] * 2 + 1 # The x dofs
            eldofs[iel,1::2] = elnods[iel,:] * 2 + 2 # The y dofs


        # Draw the mesh.
        if bDrawMesh:
            cfv.drawMesh(
                coords=coords,
                edof=eldofs,
                dofsPerNode=2,
                elType=elTypeInfo[0],
                filled=True,
                title=elTypeInfo[1])
            cfv.showAndWait()

        # Extract element coordinates
        ex, ey = cfc.coordxtr(eldofs,coords,dofs)

        # Set fixed boundary condition on left side, i.e. nodes 0-nNody
        bc = np.array(np.zeros(numNodesY*2),'i')
        idof = 1
        for i in range(numNodesY):
            idx = i*2
            bc[idx]   = idof
            bc[idx+1] = idof+1
            idof += 2

        # Assemble stiffness matrix

        K = np.zeros((ndofs,ndofs))
        R = np.zeros((ndofs,1))

        #Set the load at the right hand edge
        for i in range(numNodesY):
            R[-(i*2+2),0] = endLoadXY[0] / numNodesY
            R[-(i*2+1),0] = endLoadXY[1] / numNodesY

        for iel in range(numElements):
            if numElementNodes == 3:
                K_el, f_el = tri.tri3e(ex[iel],ey[iel],Dmat,thickness,eq)
            elif numElementNodes == 6:
                K_el, f_el = tri.tri6e(ex[iel],ey[iel],Dmat,thickness,eq)
            elif numElementNodes == 4:
                K_el, f_el = quad.quad4e(ex[iel],ey[iel],Dmat,thickness,eq)
            elif numElementNodes == 9:
                K_el, f_el = quad.quad9e(ex[iel],ey[iel],Dmat,thickness,eq)

            cfc.assem(eldofs[iel],K,K_el,R,f_el)

        r, R0 = cfc.solveq(K,R,bc)

        if K.T.all() == K.all():
            print("JA")

        nodMiddle = numNodesY//2 +1  # Mid nod on right edge
        xC = r[-(nodMiddle*2)  ,0] # 2 dofs per node, so this is the middle dof on end
        yC = r[-(nodMiddle*2)+1,0] # 2 dofs per node, so this is the middle dof on end
        print("Displacement center node right end,  x:{:12.3e}   y:{:12.3e}".format(xC, yC))

        # Sum uf reaction forces
        R0Sum = np.zeros(2,'f')
        for i in range(0,(numNodesY*2),2):
            R0Sum[0] += R0[i  ,0]
            R0Sum[1] += R0[i+1,0]

        eqTotal = eq * area* thickness #Total load for plotting purpose
        print("Total reaction force in x:{:12.3e} y:{:12.3e})".format(R0Sum[0],R0Sum[1]))

        # Draw the displacements


        if bDrawMesh:
            disp = np.array(np.zeros((numNodes,2)),'f')
            rMax = max(abs(max(r)),abs(min(r)))
            scale = 0.15 * X[1] / rMax

            for i in range( np.size(disp,0)):
                disp[i,0] = r[i*2   ,0] * scale
                disp[i,1] = r[i*2 +1,0] * scale

            cfv.drawDisplacements(displacements=disp,
                coords=coords,
                edof=eldofs,
                dofsPerNode=2,
                elType=elTypeInfo[0],
                title=elTypeInfo[1])

            cfv.showAndWait()

        dispx.append(xC)
        #dispy.append(yC)
        #elementtype.append(n)
        nodesx.append(Nx)
        #nodesy.append(Ny)


    plt.plot(nodesx,dispx)
    plt.suptitle(n)
    plt.axis([2,51,1.47e-5,1.522e-5])
    name = str(n) + str(".png")
    plt.savefig(name)
    plt.close()
    nodesx = []
    dispx=[]


