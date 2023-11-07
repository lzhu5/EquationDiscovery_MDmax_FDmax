function [A, B, uBar, cel, wavenumber, R] = generateStreamFile(H,h,T,EorS,uEorS,N)
% GENERATESTREAMFILE(H,h,T,EorS,uEorS, N)
%
% This function computes the stream function coefficients and output them
% to a file, which can be pasted into waveProperties
%
% H: Wave height [m]
% h: Water depth [m]
% T: Wave period [s]
% EorS: String. Either 'Stokes' or 'Eulerian' drift.
% uEorS: Magnitude of the Stokes or Eulerian drifts. E.g. 'Stokes' and 0
% gives a zero Stokes drift, e.g. a closed wave flume.
% N: The number of stream function coefficients. N=1 yields Airy
% theory.
%
% Modified from a stream function script by Harry Bingham, Technical
% University of Denmark.
%
% Niels G. Jacobsen, Technical University of Denmark
%

    nsteps = 500;
    g = 9.806;
    if nargin == 5
        N = 15;
    end

    % Solving for the stream function coefficients
    [eta, B, Q, c, k, R, uBar] = StreamFunctionCoefficientsPeriod(N,H,h,T,uEorS,EorS,nsteps,g);
    cel = c ;
    wavenumber = k; 

    % Solving for the surface elevation coefficients
    matrix = ones(N+1);

    for i=0:N
        matrix(i+1,1:end-1) = cos((1:N) * i * pi / N);
    end

    A = matrix \ eta';

%     % Writing the stream function coefficients in the format needed for
%     % waves2Foam
%     fid = fopen('streamFunctionCoefficient','wt');
%     fprintf(fid,' waveType\tstreamFunction;\n');
%     fprintf(fid,' N\t\t%.0f;\n',N);
%     fprintf(fid,' depth\t%f;\n',h);
%     fprintf(fid,' omega\t%f;\n',2 * pi / T);
%     fprintf(fid,' phi\t\t%f;\n',0);
%     fprintf(fid,' waveNumber\t (%f 0.0 0.0);\n',k);
%     fprintf(fid,' uBar\t%g;\n',uBar);
% 
% 
%     fprintf(fid,' A\t\tnonuniform List<scalar>\t%.0f\n (\n',N);
%     for i=1:N
%         fprintf(fid,' %g\n',A(i));
%     end
%     fprintf(fid,' );\n');
% 
% 
%     fprintf(fid,' B\t\tnonuniform List<scalar>\t%.0f\n (\n',N);
%     for i=1:N
%         fprintf(fid,' %g\n',B(i));
%     end
%     fprintf(fid,' );\n');
%     fclose(fid);

end



