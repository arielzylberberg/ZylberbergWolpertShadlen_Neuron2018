function tg = sample_tguess(tlow,thigh,N,seed)

if nargin>3 && ~isempty(seed)
    s = RandStream('twister','Seed',seed);
    use_internal_seed = 1;
else
    use_internal_seed = 0;
end
    
tlow = tlow(:);
thigh = thigh(:);

if nargin<3
    N=1;
end

if isstruct(tlow)
    thigh   = struct2array(thigh);
    tlow   = struct2array(tlow);
end

nPars = length(tlow);
if use_internal_seed 
    r = rand(s,nPars,N);
else
    r = rand(nPars,N);
end
tg = nan(nPars,N);
for i=1:N
    tg(:,i) = tlow + (thigh-tlow).*r(:,i);
end


