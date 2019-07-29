function profiles = makeProfiles(N_RUNS)

profiles = zeros(N_RUNS,4,20);

for i = 1:N_RUNS
     profiles(i,1,:) = smooth([ zeros(2,1) ; 2*rand(16,1) - 1; zeros(2,1) ],7).';
     profiles(i,2,:) = smooth([ zeros(2,1) ; 2*rand(16,1) - 1; zeros(2,1) ],7).';
     profiles(i,3,:) = smooth([ zeros(2,1) ; 2*rand(16,1) - 1; zeros(2,1) ],7).';
     profiles(i,4,:) = smooth([ zeros(2,1) ; 2*rand(16,1) - 1; zeros(2,1) ],7).';
end