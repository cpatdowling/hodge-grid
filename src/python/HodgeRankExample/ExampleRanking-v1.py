
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
import numpy
FF = QQ

#Big12Teams = ['Nebraska','Missouri','Kansas State','Kansas','Colorado','Iowa State','Oklahoma','Texas','Texas A&M','Texas Tech','Oklahoma State','Baylor']
Big12NorthTeams = ['Nebraska','Missouri','Kansas State','Kansas','Colorado','Iowa State']
Big12SouthTeams = ['Oklahoma','Texas','Texas A&M','Texas Tech','Oklahoma State','Baylor']

Big12Teams = Big12SouthTeams


Big12 = True

# load the scores
Scores = {}
Teams = set()
import re
with open('Scores','r') as fin:
  for line in fin:
    line = line.strip()
    line = [x for x in line.split(' ')[1:] if x!='']
    if '@' in line:
      line = line[:line.index('@')]
    line = ' '.join(line)
    line = [x.strip() for x in re.split('(\d+)',line) if x != '']
    Team1, Team1Score, Team2, Team2Score = line
    if Big12 and (Team1 not in Big12Teams or Team2 not in Big12Teams):
      continue
    Team1Score = int(Team1Score)
    Team2Score = int(Team2Score)
    Teams.add(Team1)
    Teams.add(Team2)
    if Team1 not in Scores.keys():
      Scores[Team1] = {}
    if Team2 not in Scores.keys():
      Scores[Team2] = {}
    Scores[Team1][Team2] = Team1Score
    Scores[Team2][Team1] = Team2Score


# remove teams that only played 1 or 2 games (this essentially recovers division 1)
TooFew = set()
for Team in Scores.keys():
  if len(Scores[Team].keys()) <=2:
    TooFew.add(Team)

for Team in TooFew:
  del Scores[Team]
  for Team1 in Scores.keys():
    if Team in Scores[Team1]: del Scores[Team1][Team]
  Teams.remove(Team)

TooFew = set()
for Team in Scores.keys():
  TooFew.add(len(Scores[Team].keys()))

# enumerate the teams
TeamToNum = {}
NumToTeam = {}
for i,Team in enumerate(Teams):
  TeamToNum[Team] = i
  NumToTeam[i] = Team

# compute pairwise score difference and construct a graph adding an edge for
# every played game
G = Graph()
ScoreDiffs = numpy.zeros((len(Teams),len(Teams)))
for Team1 in Teams:
  for Team2 in Teams:
    try:
      Diff = Scores[Team1][Team2] - Scores[Team2][Team1]
    except:
      continue
    i = TeamToNum[Team1]
    j = TeamToNum[Team2]
    ScoreDiffs[i,j] = Diff
    ScoreDiffs[j,i] = -Diff
    G.add_edge((i,j))

for Team in Teams:
  print Team, TeamToNum[Team]

## construct a graph
#G = Graph()
#Locs = numpy.vstack(numpy.where(ScoreDiffs>0)).T
#for i,j in Locs: G.add_edge((i,j))

#G.show()

print('Computing complex')
SC = G.clique_complex()
# remove all degree 3 or greater faces (not needed for now)
Faces = SC.faces()
LowDegFaces = list(Faces[0])+list(Faces[1])+list(Faces[2])
SC = SimplicialComplex(LowDegFaces)
Faces = SC.faces()
Faces = {deg:sorted([tuple(x) for x in Faces[deg]]) for deg in Faces.keys()}
MaxDeg = max(Faces.keys())

print('Computing Q')
# compute Q0, Q1, Q2
PriorRankingBig12 = {'Nebraska':172,'Missouri':124,'Kansas State':81,'Kansas':164,'Colorodo':100,'Iowa State':33,'Oklahoma':174,'Texas':174,'Texas A&M':33,'Texas Tech':89,'Oklahoma State':130,'Baylor':75}
PriorWLBig12 = {'Nebraska':10-4, 'Missouri':8-5, 'Kansas State':0, 'Iowa State':7-6, 'Colorado':3-9, 'Kansas':5-7, 'Texas':13-1, 'Oklahoma State':9-4, 'Texas Tech':9-4, 'Oklahoma':8-5, 'Texas A&M':6-7, 'Baylor':4-8}

def GetQ(TeamToNum,NumToTeam,Faces):
  Q = {}
  Q[0] = numpy.zeros((len(Faces[0]),len(Faces[0])))
  for i,u in enumerate(Faces[0]):
    for j,v in enumerate(Faces[0]):
      if j<=i:
        continue
      Q[0][i,j] = numpy.sqrt(0.01*(PriorRankingBig12[NumToTeam[u[0]]]-PriorRankingBig12[NumToTeam[v[0]]])**2 +(PriorWLBig12[NumToTeam[u[0]]]-PriorWLBig12[NumToTeam[v[0]]])**2)
      Q[0][j,i] = Q[0][i,j]
  Q[-1] = numpy.ones((1,1))
  Q[max(Faces.keys())+1] = numpy.ones((1,1))
  for dim in Faces.keys():
    if dim <=0:
      continue
    Q[dim] = numpy.zeros((len(Faces[dim]),len(Faces[dim])))
    for i,f1 in enumerate(Faces[dim]):
      for j,f2 in enumerate(Faces[dim]):
        if j <= i:
          continue
        M = numpy.zeros((len(f1),len(f1)))
        for a,u in enumerate(f1):
          for b,v in enumerate(f2):
            M[a,b] = Q[0][Faces[0].index(tuple([u])),Faces[0].index(tuple([v]))]
        Q[dim][i,j] = numpy.linalg.det(M)
        Q[dim][j,i] = Q[dim][i,j]
  return Q

Q = GetQ(TeamToNum,NumToTeam,Faces)

print('Inverting Q')
# invert Qs

Qi = {k:numpy.linalg.inv(Q[k]) for k in Q.keys()}

print('Computing differentials')
# compute d0, d1, d2
def Differentials(Faces):
  d = {}
  for k in Faces.keys():
    if k > 0:
      RangeDim = len(Faces[k-1])
      DomDim = len(Faces[k])
      Sign = 1
      d[k] = numpy.zeros((RangeDim,DomDim))
      for i,f in enumerate(Faces[k]):
        for j in range(len(f)):
          df = tuple([x for a,x in enumerate(f) if a!=j])
          if df not in Faces[k-1]:
            raise ValueError('Face not found')
          a = Faces[k-1].index(df)
          d[k][a,i]=Sign
          Sign*=-1
    d[0] = numpy.zeros((1,len(Faces[0])))
  d[max(Faces.keys())+1] = numpy.zeros((len(Faces[max(Faces.keys())]),1))
  return d

d = Differentials(Faces)

print('Computing codifferentials')
def Codifferentials(d,Q,Qi):
  return {k:numpy.dot(Qi[k],numpy.dot(d[k].T,Q[k-1])) for k in d.keys()}

delta = Codifferentials(d,Q,Qi)

Laplacians = {deg:numpy.dot(delta[deg],d[deg])+numpy.dot(d[deg+1],delta[deg+1]) for deg in Faces.keys() if deg >=0}
sys.exit(1)

print 'Computing laplacians'
Lap = Laplacians(SC,Qs,Diff,Codiff,FF)
Lap = {deg: matrix(FF,Lap[deg]) for deg in Lap.keys()}
# compute graph laplacian
L = G.laplacian_matrix()
print 'verify graph laplacian'
print Lap[0]-L


grad = Codiff[0]
curlstar = Diff[2]
Delta = Lap[1]

# compute SVD for grad
U,s,Vt = scipy.linalg.svd(grad,full_matrices=False)
V = Vt.T
Tol = 1e-12
Loc= abs(s)<=Tol
gradNullspaceBasis = V[:,Loc]
Loc = numpy.logical_not(Loc)
gradColspaceBasis = U[:,Loc]
print numpy.dot(gradNullspaceBasis.T,gradNullspaceBasis)
print numpy.dot(gradColspaceBasis.T,gradColspaceBasis)
# compute svd for curl*
U,s,Vt = scipy.linalg.svd(curlstar,full_matrices=False)
V = Vt.T
Loc= abs(s)<=Tol
curlstarNullspaceBasis = V[:,Loc]
Loc = numpy.logical_not(Loc)
curlstarColspaceBasis = U[:,Loc]
print numpy.dot(curlstarNullspaceBasis.T,curlstarNullspaceBasis)
print numpy.dot(curlstarColspaceBasis.T,curlstarColspaceBasis)
# compute svd for Delta
U,s,Vt = scipy.linalg.svd(Delta,full_matrices=False)
V = Vt.T
DeltaLoc= abs(s)<=Tol
DeltaNullspaceBasis = V[:,DeltaLoc]
DeltaLoc = numpy.logical_not(DeltaLoc)
DeltaColspaceBasis = U[:,DeltaLoc]
print numpy.dot(DeltaNullspaceBasis.T,DeltaNullspaceBasis)
print numpy.dot(DeltaColspaceBasis.T,DeltaColspaceBasis)

Basis = numpy.hstack((gradColspaceBasis,curlstarColspaceBasis,DeltaNullspaceBasis))
print numpy.dot(Basis.T,Basis)

Y = numpy.zeros(Lap[1].nrows())
for k,Edge in enumerate(Faces[1]):
  i,j = Edge
  #Y[k] = ScoreDiffs[i,j]
  Y[k] = numpy.sign(ScoreDiffs[i,j]) # using just the sign reproduces http://www.ams.org/samplings/feature-column/fc-2012-12

grads = numpy.dot(gradColspaceBasis,numpy.dot(Y.T,gradColspaceBasis))
curlstarPhi = numpy.dot(curlstarColspaceBasis,numpy.dot(Y.T,curlstarColspaceBasis))
h = numpy.dot(DeltaNullspaceBasis,numpy.dot(Y.T,DeltaNullspaceBasis))

print 'edge\tgrads\tcurlstarPhi\th'
for k, Edge in enumerate(Faces[1]):
  i,j = Edge
  print NumToTeam[i],NumToTeam[j],'\t',grads[k],'\t',curlstarPhi[k],'\t',h[k]

# compute s and phi using pseudo-inverse
s = numpy.dot(numpy.linalg.pinv(grad),grads)
Phi = numpy.dot(numpy.linalg.pinv(curlstar),curlstarPhi)

# get ker(L) basis for C_{0} and C_{1}
# compute svd for L0
U,Sigma,Vt = scipy.linalg.svd(Lap[0],full_matrices=False)
V = Vt.T
L0Loc= abs(Sigma)<=Tol
L0NullspaceBasis = V[:,L0Loc]
L0Loc = numpy.logical_not(L0Loc)
L0ColspaceBasis = U[:,L0Loc]
print numpy.dot(L0NullspaceBasis.T,L0NullspaceBasis)
print numpy.dot(L0ColspaceBasis.T,L0ColspaceBasis)
# compute svd for L0
U,Sigma,Vt = scipy.linalg.svd(Lap[2],full_matrices=False)
V = Vt.T
L2Loc= abs(Sigma)<=Tol
L2NullspaceBasis = V[:,L2Loc]
L2Loc = numpy.logical_not(L2Loc)
L2ColspaceBasis = U[:,L2Loc]
print numpy.dot(L2NullspaceBasis.T,L2NullspaceBasis)
print numpy.dot(L2ColspaceBasis.T,L2ColspaceBasis)

# project out the harmonic piece to get the distinguished harmonic rep
sharmonic = numpy.dot(L0NullspaceBasis,numpy.dot(s.T,L0NullspaceBasis))
Phiharmonic = numpy.dot(L2NullspaceBasis,numpy.dot(Phi.T,L2NullspaceBasis))
sDistinguished = s-sharmonic
PhiDistinguished = Phi-Phiharmonic

## compute function
f = {deg: numpy.zeros(Lap[deg].nrows()) for deg in Lap.keys()}
for k,Edge in enumerate(Faces[1]):
  i,j = Edge
  f[1][k] = ScoreDiffs[i,j]
print 'compute decomposition'
Imd,Imdelta,KerL = HodgeDecomposition(f,Diff,Codiff,Lap,Qs,Tol=1e-12)
#df1,f1 = Imd[1]
#deltaf2,f2 = Imdelta[1]
#f3 = KerL[1]
#
#h1 = KerL[1]
grads1,s1 = Imdelta[1]
sDistinguished = s1
#curlstarPhi1, Phi1 = Imd[1]
#
#Ranking = numpy.argsort(s1)
#print 'ranking'
#for k in Ranking:
#  print NumToTeam[k]



Big12Teams = ['Nebraska','Missouri','Kansas State','Kansas','Colorado','Iowa State','Oklahoma','Texas','Texas A&M','Texas Tech','Oklahoma State','Baylor']
Big12NorthTeams = ['Nebraska','Missouri','Kansas State','Kansas','Colorado','Iowa State']
Big12NorthTeams = ['Nebraska','Missouri','Kansas State','Kansas','Colorado','Iowa State']
Big12SouthTeams = ['Oklahoma','Texas','Texas A&M','Texas Tech','Oklahoma State','Baylor']
sBig12 = [sDistinguished[TeamToNum[x]] for x in Big12Teams]
sBig12N = [sDistinguished[TeamToNum[x]] for x in Big12NorthTeams]
sBig12S = [sDistinguished[TeamToNum[x]] for x in Big12SouthTeams]
Big12Ranking = numpy.argsort(sBig12)
Big12NorthRanking = numpy.argsort(sBig12N)
Big12SouthRanking = numpy.argsort(sBig12S)
Big12TeamRanking = [Big12Teams[x] for x in Big12Ranking]
Big12NorthTeamRanking = [Big12NorthTeams[x] for x in Big12NorthRanking]
Big12SouthTeamRanking = [Big12SouthTeams[x] for x in Big12SouthRanking]
TotalRanking = numpy.argsort(sDistinguished)
RankingScores = sDistinguished.copy()#[TotalRanking]
RankingScores -= RankingScores.min()
RankingScores /= RankingScores.max()
RankingScores *= 100
RankingScores = 100-RankingScores
Top5 = TotalRanking[:5]
Bottom5 = TotalRanking[-5:]
print 'Top 5:'
for k in Top5:
  print NumToTeam[k], RankingScores[k]
print
print 'Bottom 5:'
for k in Bottom5:
  print NumToTeam[k]

print 'Big 12', Big12TeamRanking
print 'Big 12 North', Big12NorthTeamRanking
print 'Big 12 South', Big12SouthTeamRanking

#
#
#
#
#
#
#
#
#
#

## compute function
#f = {deg: numpy.zeros(Lap[deg].nrows()) for deg in Lap.keys()}
#for k,Edge in enumerate(Faces[1]):
#  i,j = Edge
#  f[1][k] = ScoreDiffs[i,j]
#print 'compute decomposition'
#Imd,Imdelta,KerL = HodgeDecomposition(f,Diff,Codiff,Lap,Qs,Tol=1e-12)
#df1,f1 = Imd[1]
#deltaf2,f2 = Imdelta[1]
#f3 = KerL[1]
#
#h1 = KerL[1]
#grads1,s1 = Imdelta[1]
#curlstarPhi1, Phi1 = Imd[1]
#
#Ranking = numpy.argsort(s1)
#print 'ranking'
#for k in Ranking:
#  print NumToTeam[k]
