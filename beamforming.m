
% This example introduces the basic concept of hybrid beamforming 
% and shows how to split the precoding and combining weights using orthogonal matching pursuit algorithm. 
% It shows that hybrid beamforming can closely match the performance offered by optimal digital weights.

% This program simulates a 64 x 16 MIMO hybrid beamforming system, 
% with a 64-element square array with 4 RF chains on the transmitter side 
% and a 16-element square array with 4 RF chains on the receiver side.
% Each antenna element can be connected to one or more TR modules.

Nt = 64;    % the number of antenna elements on the transmit side
NtRF = 4;   % the number of TR switches on the transmit side

Nr = 16;    % the number of antenna elements on the receive side
NrRF = 4;   % the number of TR switches on the receive side

% it is assumed that each antenna is connected to all RF chains. Thus, each antenna is connected to 4 phase shifters. 
% Such an array can be modeled by partitioning the array aperture into 4 completely connected subarrays.

rng(4096);
c = 3e8;    % Speed of Light
fc = 28e9;  % Frequency
lambda = c/fc;  % Wavelength

% Create transmit uniform rectangular array partionned into subarrays whose elements are spaced half a wavelength apart.
txarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nt) sqrt(Nt)],lambda/2),...
    'SubarraySelection',ones(NtRF,Nt),'SubarraySteering','Custom'); % subarrays are steered by setting independent weights for all elements in each subarray.
% Create receive uniform rectangular array partionned into subarrays whose elements are spaced half a wavelength apart.
rxarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nr) sqrt(Nr)],lambda/2),...
    'SubarraySelection',ones(NrRF,Nr),'SubarraySteering','Custom');



% The path gain for each scatterer is obtained from a complex circular symmetric Gaussian distribution.
Ncl = 6;    % assume a scattering environment with 6 scattering clusters randomly distributed in space. 
Nray = 8;   % Within each cluster, there are 8 closely located scatterers. 
Nscatter = Nray*Ncl;    % a total of 48 scatterers
angspread = 5;  % Within each cluster, there are closely located scatterers with a angle spread of 5 degrees. 

% compute randomly placed scatterer clusters
txclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30];
rxclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30];
txang = zeros(2,Nscatter);
rxang = zeros(2,Nscatter);
% compute the rays within each cluster
for m = 1:Ncl
    txang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+txclang(:,m);
    rxang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+rxclang(:,m);
end

g = (randn(1,Nscatter)+1i*randn(1,Nscatter))/sqrt(Nscatter);

% The channel matrix can be formed as : 
txpos = getElementPosition(txarray)/lambda; % Get transmitter position
rxpos = getElementPosition(rxarray)/lambda; % Get receiver position
H = scatteringchanmtx(txpos,rxpos,txang,rxang,g);

% Hybrid Weights Computation

% Assuming the channel is known, the unconstrained optimal precoding weights can be obtained by diagonalizing 
% The channel matrix and extracting the first NtRF dominating modes. 

F = diagbfweights(H);
F = F(1:NtRF,:);

% The transmit beam pattern can be plotted as:
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',F','PropagationSpeed',c);
% Even in a multipath environment, there are limited number of dominant directions.

% The hybrid weights can be computed as:
At = steervec(txpos,txang);
Ar = steervec(rxpos,rxang);

Ns = NtRF;
[Fbb,Frf] = omphybweights(H,Ns,NtRF,At);

% The beam pattern of the hybrid weights:
pattern(txarray,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',Frf'*Fbb','PropagationSpeed',c);

% Spectral Efficiency Comparison
% Compares the spectral efficiency achieved using the optimal weights with that of the proposed hybrid beamforming weights. 

%The resulting spectral efficiency curve is obtained from 50 Monte-Carlo trials for each SNR.
snr_param = -40:5:0;
Nsnr = numel(snr_param);
Ns_param = [1 2];
NNs = numel(Ns_param);

NtRF = 4;
NrRF = 4;

Ropt = zeros(Nsnr,NNs);
Rhyb = zeros(Nsnr,NNs);
Niter = 50;

% The transmit antenna array is assumed to be at a base station,
% with a focused beamwidth of 60 degrees in azimuth and 20 degrees in elevation. 
for m = 1:Nsnr
    snr = db2pow(snr_param(m));
    for n = 1:Niter
        % Channel realization
        txang = [rand(1,Nscatter)*60-30;rand(1,Nscatter)*20-10];
        rxang = [rand(1,Nscatter)*180-90;rand(1,Nscatter)*90-45];
        At = steervec(txpos,txang);
        Ar = steervec(rxpos,rxang);
        g = (randn(1,Nscatter)+1i*randn(1,Nscatter))/sqrt(Nscatter);
        H = scatteringchanmtx(txpos,rxpos,txang,rxang,g);
        
        for k = 1:NNs
            Ns = Ns_param(k);
            % Compute optimal weights and its spectral efficiency
            [Fopt,Wopt] = helperOptimalHybridWeights(H,Ns,1/snr);
            Ropt(m,k) = Ropt(m,k)+helperComputeSpectralEfficiency(H,Fopt,Wopt,Ns,snr);

            % Compute hybrid weights and its spectral efficiency
            [Fbb,Frf,Wbb,Wrf] = omphybweights(H,Ns,NtRF,At,NrRF,Ar,1/snr);
            Rhyb(m,k) = Rhyb(m,k)+helperComputeSpectralEfficiency(H,Fbb*Frf,Wrf*Wbb,Ns,snr);
        end
    end
end
Ropt = Ropt/Niter;
Rhyb = Rhyb/Niter;

% the hybrid beamforming can perform close to what optimal weights can offer using less hardware.
plot(snr_param,Ropt(:,1),'--sr',...
    snr_param,Ropt(:,2),'--b',...
    snr_param,Rhyb(:,1),'-sr',...
    snr_param,Rhyb(:,2),'-b');
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz');
legend('Ns=1 optimal','Ns=2 optimal','Ns=1 hybrid', 'Ns=2 hybrid',...
    'Location','best');
grid on;