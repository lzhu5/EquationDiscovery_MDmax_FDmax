% clc


H = 0.18; 
h = 0.4 ;
T = 2.3 ;

grav = 9.806 ;

EorS = 'Eulerian' ;
uEorS = 0.0 ;
N = 15 ; 

[A, B, uBar, cel, wavenumber] = generateStreamFile(H,h,T,EorS,uEorS,N) ; 
k = wavenumber ; 
L = 2*pi/k 

return

fid_A = fopen ('Acoeff.txt', 'w') ; 
fid_B = fopen ('Bcoeff.txt', 'w') ; 
fid_uBar = fopen ('uBar.txt', 'w') ; 
fid_cel = fopen ('cel.txt', 'w') ; 
fid_wvnum = fopen ('wavenumber.txt', 'w') ; 

fprintf (fid_A, '%20.10f ', A) ; 
fprintf (fid_B, '%20.10f ', B) ; 
fprintf (fid_uBar, '%20.10f ', uBar) ; 
fprintf (fid_cel, '%20.10f ', cel) ; 
fprintf (fid_wvnum, '%20.10f ', wavenumber) ; 

fclose (fid_A) ;
fclose (fid_B) ;
fclose (fid_uBar) ;
fclose (fid_cel) ;
fclose (fid_wvnum) ;
 
