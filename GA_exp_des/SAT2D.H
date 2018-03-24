C
C---------------------------------------------------------------------72
C
C Dimensioning Parameters for Sat3d Code
C --------------------------------------
C   NMAX   - maximum N
C            N     = # of nodes
C   NTRMAX - maximum NTRI  
C            NTRI  = # of triangles
C   NPMAX  - maximum NP
C            NDIR  = # of Dirichlet nodes
C   NQMAX  - maximum NQ    
C            NQ    = # of Neumann nodes in 2-d mesh
C   MAXZON - maximum NZONE 
C            NZONE = # of material types in the porous medium
C   MAXTRM - maximum NTERM 
C            NTERM = # of nonzero elements in system matrices
C   NRMAX  - maximum NR    
C            NR    = # of nodes selected for partial output
C   MAXPRT - maximum NPRT  
C            NPRT  = # of time values for detailed output
C   N1MAX  - maximum N1    
C            N1    = maximum # of element connections to a node
C   NTPMAX - N1MAX*NMAX
C
C---------------------------------------------------------------------72
C
      INTEGER   NMAX,NTRMAX
      INTEGER   NPMAX,NQMAX
      INTEGER   MAXZON,MAXTRM
      INTEGER   NRMAX,MAXPRT,N1MAX,NTPMAX
      PARAMETER (NMAX=29197,NTRMAX=57888)
      PARAMETER (NPMAX=504,NQMAX=25)
      PARAMETER (NRMAX=25,MAXPRT=200,N1MAX=20,NTPMAX=N1MAX*NMAX)
      PARAMETER (MAXZON=15,MAXTRM=N1MAX*NMAX)
