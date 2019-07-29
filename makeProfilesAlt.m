function profiles = makeProfilesAlt(N_RUNS)

profiles = zeros(N_RUNS,4,100);

for i = 1:N_RUNS
     profiles(i,1,:) = zeros(100,1);
     profiles(i,2,:) = smooth([ zeros(25,1) ; 2*rand(50,1) - 1; zeros(25,1) ],7).';
     profiles(i,3,:) = zeros(100,1);
     profiles(i,4,:) = zeros(100,1);
end