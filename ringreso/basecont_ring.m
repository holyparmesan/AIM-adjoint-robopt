function baseconts = basecont_ring(spacing)

% ring geometry
R_in = 4.1e-6;
R_out = 4.6e-6;
y = 5e-6;
%y = 0; % can't do this - if you try to factor it in, you're going to get lost finding
%reference positions for adjointifying

% the understanding is t_l goes from 0 to 1.
N = 2*pi*R_out / spacing;
t_l = linspace(0,1,N);

xs_out = R_out * sin(2*pi*t_l);
ys_out = y - R_out * cos(2*pi*t_l);
xs_in = R_in * sin(2*pi*t_l);
ys_in = y - R_in * cos(2*pi*t_l);

baseconts = cell(2,1);
baseconts{1} = [xs_out ; ys_out].';
baseconts{2} = [xs_in ; ys_in].';

end