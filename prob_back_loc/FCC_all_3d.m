function [ env, tlag1, pair1 ] = FCC_all_3d( Nsta, sig, lag )
% =========================================================================
% Function: FCC_all.m
%
% Calculate first cross-correlation (FCC) between all pairs, and return
% central part of FCC envelopes based on defined lag.
%
% For real-data example
% =========================================================================
t1 = tic;

[ pair1, Npair1 ] = cons_pair1( Nsta );

[n,m] = size(sig);
tlag = -(m-1):(m-1);
env =  zeros(Npair1,2*lag+1);

for ii = 1:size(pair1,1)
    %sig_c = xcorr(sig(pair1(ii,1),:),sig(pair1(ii,2),:));
    sig_c = xcorr(sig(pair1(ii,1),:),sig(pair1(ii,2),:),'coeff');
    sig_env2 = abs(hilbert(sig_c));
    env(ii,:) = sig_env2(1,find(tlag==-lag):find(tlag==lag));
    clear sig_c sig_env2
end

fprintf('Time used for computing FCCs and envelopes = %f s\n',toc(t1));
tlag1 = -lag:lag;
