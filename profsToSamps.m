function samps = profsToSamps(profiles, multipliers)

N_MULT = length(multipliers);
N_PROF = size(profiles,1);

sz = size(profiles);
sz(1) = N_MULT*N_PROF;
samps = zeros(sz);
% we pretend we don't know that profiles is 3-dimensional

for i = 1:N_PROF
  for j = 1:N_MULT % and then require it here YOLO
    samps(N_MULT*(i-1) + j,:,:) = profiles(i,:,:) * multipliers(j);
  end
end

end