#! /usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Paul Bruillard'
__copyright__ = 'Copyright 2014, Paul Bruillard'
__license__ = 'MIT License\n\
\n\
Copyright (C) 2014 Paul Bruillard\n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n\
\n\
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\
\n\
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'
__version__ = '0.0'
__maintainer__ = 'Paul Bruillard'
__email__ = 'Paul.Bruillard@pnnl.gov'
__status__ = 'Development'
__date__ = ''
__credits__ = ''

def ComputeAdjoint(A,Qdom,QIm,FF):
  '''
  Returns the adjoint of A with respect to the inner products on the domain and
  range given by Qdom and QIm respectively.
  '''
  # Recall that Aad is defined by xdagger*Aaddagger*Qdom*y(Aad x,y) = (x,A y) =
  # xdagger*Qim*A*y
  # this is true for all x and y so we have
  # Aaddagger = Qim*A*Qdom.inv
  A = matrix(FF,A)
  Qdom = matrix(FF,Qdom)
  QIm = matrix(FF,QIm)
  return (QIm*A*(Qdom.inverse())).T.conjugate()

import scipy.linalg
def perm_parity(lst):
  '''\
  Given a permutation of the digits 0..N in order as a list, 
  returns its parity (or sign): +1 for even parity; -1 for odd.
  '''
  parity = 1
  SortedList = sorted(lst)
  for i in range(0,len(lst)-1):
    if lst[i] != SortedList[i]:
      parity *= -1
      mn = min(range(i,len(lst)), key=lst.__getitem__)
      lst[i],lst[mn] = lst[mn],lst[i]
  return parity  

# A function to produce differentials for the simplicial complex
#
# SAGE does this, but I've been unable to get the basis in terms of simplicies
# back so there is no easy way to apply Q.
def DifferentialElement(FaceSimplex,LowerBasis):
  Face = list(FaceSimplex)
  LowerBasisList = list(LowerBasis)
  d = [0 for _ in LowerBasis]
  Sign = 1
  if len(LowerBasis)==1 and len(list(LowerBasis[0]))==0 or len(Face)==0:
    return d[:]

  for k in xrange(len(Face)):
    i = LowerBasisList.index(Simplex(Face[:k]+Face[k+1:]))
    d[i] += Sign
    Sign *= -1
  return d[:]

def DifferentialMap(LargerBasis,SmallerBasis,FF):
  return matrix(FF,[DifferentialElement(FaceSimplex,SmallerBasis) for FaceSimplex in LargerBasis]).T

def Differentials(Complex,FF):
  Faces = Complex.faces()
  Faces = {deg:sorted(Faces[deg]) for deg in Faces.keys()}
  MaxDeg = max(Faces.keys())
  Diffs = {deg: DifferentialMap(Faces[deg],Faces[deg-1],FF) for deg in xrange(MaxDeg+1)}
  Diffs[MaxDeg+1] = matrix(FF,len(Faces[MaxDeg]),1,0)
  Diffs[MaxDeg+2] = matrix(FF,1,1,0)
  return Diffs

# A function to compute 
#  Q((a_{1},...,a_{k}),(b_{1},...,b_{k}))=det(Q[a_{i},b_{j}])/k!
# this is used to induce a vertex form up to higher degree
def BSQ(InitialSimplex,TerminalSimplex,Q0,QLast,Faces0,FacesLast,FF):
  i = FacesLast.index(InitialSimplex[1:])
  i0 = Faces0.index([InitialSimplex[0]])
  Sign = 1
  Val = 0
  for k in xrange(len(TerminalSimplex)):
    if k == 0:
      j = FacesLast.index(TerminalSimplex[1:])
      SubTerm = TerminalSimplex[1:]
    else:
      SubTerm = TerminalSimplex[:k]+TerminalSimplex[k+1:]
      j = FacesLast.index(TerminalSimplex[:k]+TerminalSimplex[k+1:])
    j0 = Faces0.index([TerminalSimplex[k]])
#    print '\t', TerminalSimplex[k],SubTerm, j0, j
    Val+=Q0[i0,j0]*QLast[i,j]*Sign
    Sign*=-1
  return Val

def InduceQ(Q,Complex,FF):
  # induce up the inner product if needed
  Qs = {}
  Faces = Complex.faces()
  Faces = {deg: [list(x) for x in sorted(Faces[deg])] for deg in Faces.keys()}
  MaxDeg = max(Faces.keys())
  Qs[0] = Q
  for dim in Faces.keys():
    if dim <= 0:
      continue
    if dim > MaxDeg:
      break
    Qs[dim] = matrix(FF,len(Faces[dim]),len(Faces[dim]),0)
    for InitialCounter in xrange(len(Faces[dim])):
      InitialFace = Faces[dim][InitialCounter]
      Qs[dim][InitialCounter,InitialCounter] = BSQ(InitialFace,InitialFace,Q,Qs[dim-1],Faces[0],Faces[dim-1],FF)
      for FinalCounter in xrange(InitialCounter+1,len(Faces[dim])):
        TerminalFace = Faces[dim][FinalCounter]
        Qs[dim][InitialCounter,FinalCounter] = BSQ(InitialFace,TerminalFace,Q,Qs[dim-1],Faces[0],Faces[dim-1],FF)
        Qs[dim][FinalCounter,InitialCounter] = Qs[dim][InitialCounter,FinalCounter]
    Qs[dim-1]/=factorial(dim)
  Qs[MaxDeg]/=factorial(MaxDeg+1)
  Qs[-1] = matrix(FF,1,1,0)
  Qs[-2] = matrix(FF,1,1,0)
  Qs[MaxDeg+1] = matrix(FF,1,1,0)
  Qs[MaxDeg+2] = matrix(FF,1,1,0)
  return Qs

def GetQs(Q,Complex,FF):
  '''
  Computes the Qs and QsInv for the Complex.

  If Q='Id', then all Qs are the identity.
  If Q='r', then Q is randomized and induced up.
  If Q is a matrix, then it is induced up.

  WARNING: Each Q is inverted. This can take a significant amount of time.
  '''
  Faces = Complex.faces()
  Faces = {deg:sorted(Faces[deg]) for deg in Faces.keys()}
  print 'Computing Q'
  if Q!='Id':
    if Q == 'r':
      # select a nondegenerate hermitian form (on vertices)
      # first we generate a random hermitian matrix which is nondegenerate
      Q = numpy.random.rand(NumNodes,NumNodes)
      Q += Q.T.conj()
      Q = numpy.dot(Q,Q.T.conj())
      while numpy.linalg.det(Q)==0:
        Q = numpy.random.rand(NumNodes,NumNodes)
        Q += Q.T.conj()
        # multiply by hermitian conjugate to get positive semidefinite
        Q = numpy.dot(Q,Q.T.cong())
      Q = matrix(FF,Q)
    print 'Inducing Q'
    Qs = InduceQ(Q,CliqueComplex,FF)
    print 'Inverting Q'
    for deg in Qs.keys():
      print deg,Qs[deg].nrows(), Qs[deg].base_ring()
    QsInv = {deg: Qs[deg].inverse() if Qs[deg].det()!=0 else Qs[deg] for deg in Qs.keys()}
  else:
    for deg in Faces.keys():
      print len(Faces[deg])
    Qs = {deg:matrix.identity(len(Faces[deg])) for deg in Faces.keys()}
    QsInv = {deg:matrix.identity(len(Faces[deg])) for deg in Faces.keys()}
    Qs[-1] = matrix(FF,1,1,0)
    Qs[-2] = matrix(FF,1,1,0)
    Qs[MaxDeg+1] = matrix(FF,1,1,0)
    Qs[MaxDeg+2] = matrix(FF,1,1,0)
    QsInv[-1] = matrix(FF,1,1,0)
    QsInv[-2] = matrix(FF,1,1,0)
    QsInv[MaxDeg+1] = matrix(FF,1,1,0)
    QsInv[MaxDeg+2] = matrix(FF,1,1,0)
  return Qs, QsInv

def Codifferentials(Diff,Qs,QsInv,FF):
  MaxDeg = max(Qs.keys())-2
  Codiffs = {deg: QsInv[deg+1]*(Diff[deg+1].H)*Qs[deg] for deg in xrange(MaxDeg)}
  Codiffs[-1] = Diff[0].T.conjugate()*0
  Codiffs[-2] = matrix(FF,1,1,0)
  Codiffs[MaxDeg] = matrix(FF,Diff[MaxDeg+1].ncols(),Qs[MaxDeg].ncols(),0)
  Codiffs[MaxDeg+1] = matrix(FF,1,1,0)
  return Codiffs

# A function to compute Laplacians for a simplicial complex to all orders and
# using a specified hermitian form.
def Laplacians(Complex,Qs,d,delta,FF):
  '''
  Complex should be a SAGE simplicial complex and Qs should either be a
  dictionary specifying a hermitian form on each degree or a matrix specifying a
  hermtian inner product on the lowest degree. In the latter case the inner
  product will be induced up.

  FF is the base ring that the matrices are constructed over. This should be a
  field and it should match the base ring that the complex will be defined over.
  It should be possible to compute
  Q((a_{1},...,a_{k}),(b_{1},...,b_{k}))=det(Q[a_{i},a_{j}])/k! for all k-faces
  over this field.
  '''
  # get the simplicies in the complex
  Faces = Complex.faces()

  L = {}
  for dim in Faces.keys():
    L[dim] = matrix(FF,len(Faces[dim]),len(Faces[dim]),0)
    if dim+1 in d.keys() and dim in delta.keys():
      L[dim] += d[dim+1]*delta[dim]
    if dim-1 in delta.keys() and dim in d.keys():
      L[dim]+=delta[dim-1]*d[dim]
  return L

def GS(Min,Q):
  '''performs gram-schmidt orthogonalization on the columns of M under the inner
  product defined by Q.'''
  M = Min.copy()
  for k in xrange(M.shape[1]):
    Norm = numpy.sqrt(numpy.dot(M[:,k].T,numpy.dot(Q,M[:,k])))
    if Norm == 0:
      M[:,k] *= 0
      continue
    else:
      M[:,k]/=Norm
    for l in xrange(k+1,M.shape[1]):
      M[:,l]-=numpy.dot(M[:,l].T,numpy.dot(Q,M[:,k]))*M[:,k]
  # remove columns of zeros
  M = M.compress(numpy.logical_not(numpy.all(M==0,axis=0)),axis=1)
  return M

def KernelImg(M,QNull,QCol,Tol=1e-6):
  # computes a basis for the image (column space) and kernel (nullspace) of M
  # subject to the tolerance :Tol:. These bases must be orthonormal under the
  # Q-inner product
  U,s,Vt = scipy.linalg.svd(M,full_matrices=False)
  V = Vt.T
  Loc = abs(s)<=Tol
  NullspaceBasis = V[:,Loc]
  Loc = numpy.logical_not(Loc)
  ColspaceBasis = U[:,Loc]

  # orthonormalize each under the Q-inner product
  if NullspaceBasis.size > 0 and numpy.abs(numpy.dot(NullspaceBasis.T,numpy.dot(QNull,NullspaceBasis))-numpy.eye(NullspaceBasis.shape[1])).max() > Tol:
    NullspaceBasis = GS(NullspaceBasis,QNull)
  if ColspaceBasis.size > 0 and numpy.abs(numpy.dot(ColspaceBasis.T,numpy.dot(QCol,ColspaceBasis))-numpy.eye(ColspaceBasis.shape[1])).max() > Tol:
    ColspaceBasis = GS(ColspaceBasis,QCol)
  # verify orthogonality
  if NullspaceBasis.size > 0 and numpy.abs(numpy.dot(NullspaceBasis.T,numpy.dot(QNull,NullspaceBasis))-numpy.eye(NullspaceBasis.shape[1])).max() > Tol:
    sys.stderr.write('Nullspace error\n')
    sys.stderr.flush()
    sys.exit(1)
  if ColspaceBasis.size > 0 and numpy.abs(numpy.dot(ColspaceBasis.T,numpy.dot(QCol,ColspaceBasis))-numpy.eye(ColspaceBasis.shape[1])).max() > Tol:
    sys.stderr.write('Colspace error\n')
    sys.stderr.flush()
    sys.exit(1)
  return NullspaceBasis,ColspaceBasis

def HodgeDecomposition(f,Diff,Codiff,Lap,Qs,Tol=1e-6):
  # compute the nullspace and column space of Diff, Codiff, and Lap
  LapNullSpaces = {}
  DiffNullSpaces = {}
  CodiffNullSpaces = {}
  LapColSpaces = {}
  DiffColSpaces = {}
  CodiffColSpaces = {}
  for deg in f.keys():
    print 'lap spaces in deg', deg
    LapNullSpaces[deg],LapColSpaces[deg] = KernelImg(Lap[deg],Qs[deg],Qs[deg],Tol)
    if deg in Diff.keys():
      print 'diff spaces in deg', deg
      DiffNullSpaces[deg],DiffColSpaces[deg] = KernelImg(Diff[deg],Qs[deg],Qs[deg-1],Tol)
    elif deg not in Diff.keys():
      DiffNullSpaces[deg] = numpy.array([[0]])
      DiffColSpaces[deg] = numpy.zeros((f[deg-1].shape[0] if deg-1 in f.keys() else 1,1))
    if deg in Codiff.keys():
      print 'codiff spaces in deg', deg
      CodiffNullSpaces[deg],CodiffColSpaces[deg] = KernelImg(Codiff[deg],Qs[deg],Qs[deg+1],Tol)
    elif deg not in Codiff.keys():
      CodiffNullSpaces[deg] = numpy.array([[0]])
      CodiffColSpaces[deg] = numpy.zeros((f[deg+1].shape[0] if deg+1 in f.keys() else 1,1))
    # verify that the differential and codifferential annihilate the laplacian
    # nullspace
    if deg in Diff.keys() and LapNullSpaces[deg].size > 0 and numpy.abs(numpy.dot(Diff[deg],LapNullSpaces[deg])).max() > Tol:
      sys.stderr.write('Diff does not annihilate ker(Lap) in degree '+str(deg)+'\n')
      sys.stderr.flush()
      sys.exit(1)
    if deg in Codiff.keys() and LapNullSpaces[deg].size > 0 and numpy.abs(numpy.dot(Codiff[deg],LapNullSpaces[deg])).max() > Tol:
      sys.stderr.write('Codiff does not annihilate ker(Lap) in degree '+str(deg)+'\n')
      sys.stderr.flush()
      sys.exit(1)
  MaxDeg = max(f.keys())
  DiffColSpaces[MaxDeg+1] = numpy.zeros((f[MaxDeg].shape[0],0))
  DiffNullSpaces[MaxDeg+1] = numpy.zeros((1,0))
  CodiffColSpaces[-2] = numpy.zeros((1,0))
  CodiffNullSpaces[-2] = numpy.zeros((1,0))

  # in each degree verify that the im, im, ker basis is orthonormal with respect
  # to Q
  for deg in f.keys():
    dIm = DiffColSpaces[deg+1]
    deltaIm = CodiffColSpaces[deg-1]
    Lker = LapNullSpaces[deg]
    print dIm.shape, deltaIm.shape, Lker.shape
    Basis = numpy.hstack((dIm,deltaIm,Lker))
    #Basis = numpy.hstack((deltaIm,Lker))
    NormdIm = numpy.dot(dIm.T,numpy.dot(Qs[deg],dIm))
    NormdIm -= numpy.eye(NormdIm.shape[0])
    if NormdIm.size >0:
      NormdImMax = numpy.abs(NormdIm).max()
    else:
      NormdImMax = 0
    NormdeltaIm = numpy.dot(deltaIm.T,numpy.dot(Qs[deg],deltaIm))
    NormdeltaIm -= numpy.eye(NormdeltaIm.shape[0])
    if NormdeltaIm.size > 0:
      NormdeltaImMax = numpy.abs(NormdeltaIm).max()
    else:
      NormdeltaImMax = 0
    NormLKer = numpy.dot(Lker.T,numpy.dot(Qs[deg],Lker))
    NormLKer -= numpy.eye(NormLKer.shape[0])
    NormLKerMax = numpy.abs(NormLKer).max() if NormLKer.size>0 else 0
    Norm = numpy.dot(Basis.T,numpy.dot(Qs[deg],Basis))
    Norm -= numpy.eye(Norm.shape[0])
    NormMax = numpy.abs(Norm).max() if Norm.size>0 else 0
    print deg
    print '\t d', NormdImMax
    print '\t d*', NormdeltaImMax
    print '\t L', NormLKerMax
    print '\t', NormMax


  Imd = {}
  Imdelta = {}
  KerL = {}
  for deg in f.keys():
    # get kernel piece, and two image pieces
    NullspaceCoeff = numpy.dot(f[deg],numpy.dot(Qs[deg],LapNullSpaces[deg]))
    CodiffImgCoeff = numpy.dot(f[deg],numpy.dot(Qs[deg],CodiffColSpaces[deg-1]))
    DiffImgCoeff = numpy.dot(f[deg],numpy.dot(Qs[deg],DiffColSpaces[deg+1]))
    # verify reconstruction
    f3 = numpy.dot(LapNullSpaces[deg],NullspaceCoeff)
    deltaf2 = numpy.dot(CodiffColSpaces[deg-1],CodiffImgCoeff)
    df1 = numpy.dot(DiffColSpaces[deg+1],DiffImgCoeff)
    Recon = df1+deltaf2+f3

    if numpy.abs(Recon-f[deg]).max() > Tol:
      sys.stderr.write('Error recon 1 in degree '+str(deg)+'\n')
      sys.stderr.flush()
      sys.exit(1)
    # compute f1 and f2. Note that f1 is only defined up to an element of ker(d)
    # and f2 is only defined up to an element of ker(delta). To avoid any
    # ambiguity we will project out the kernel
    f1 = numpy.dot(numpy.linalg.pinv(Diff[deg+1]),df1)
    f1ker = numpy.dot(numpy.dot(f1,numpy.dot(Qs[deg+1],DiffNullSpaces[deg+1])),DiffNullSpaces[deg+1].T)
    f1 = f1-f1ker
    f2 = numpy.dot(numpy.linalg.pinv(Codiff[deg-1]),deltaf2)
    f2ker = numpy.dot(numpy.dot(f2,numpy.dot(Qs[deg-1],CodiffNullSpaces[deg-1])),CodiffNullSpaces[deg-1].T)
    f2 = f2-f2ker
    # verify reconstruction
    Recon = numpy.dot(Diff[deg+1],f1)+numpy.dot(Codiff[deg-1],f2)+f3
    if numpy.abs(Recon-f[deg]).max() > Tol:
      sys.stderr.write('Error recon 2\n')
      sys.stderr.flush()
      sys.exit(1)
    KerL[deg] = f3
    Imd[deg] = (df1,f1)
    Imdelta[deg] = (deltaf2,f2)
  return Imd,Imdelta,KerL

