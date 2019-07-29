function conts = yb_contours(spacing)

% returns 4 contours: top out, top in, bottom in, bottom out

% We have decided that the flat sections after the bends extend for 1um,
% from 12.5 to 13.5, before we stop caring about the effects of roughness
% and that before the input, it extends from -1.5um to 0um

N_BENDI = floor(10e-6 / spacing);
N_BENDO = floor(10e-6 / spacing);
N_CAV = floor(3 * 1.5e-6 / spacing);
N_INPUT = floor(1.5e-6 / spacing);
N_OUTPUT = floor(1e-6 / spacing);

[ys_inbend, xs_inbend] = yb_in(linspace(0,1,N_BENDI));
[ys_out,xs_out] = yb_out(linspace(0,1,N_BENDO));
[ys_cav,xs_cav] = yb_cav(linspace(0,1,N_CAV));

ys_base = ones(1,N_INPUT) * 0.25;
xs_base = linspace(-1.5,0,N_INPUT);

xs_ext = linspace(12.5,13.5,N_OUTPUT);
ys_inext = ones(1,N_OUTPUT) * 2.5;
ys_outext = ones(1,N_OUTPUT) * 3;

top_out = [xs_base xs_cav(2:end) xs_out(2:end) xs_ext(2:end);ys_base ys_cav(2:end) ys_out(2:end) ys_outext(2:end)].'*1e-6;
bot_out = [top_out(:,1) -1*top_out(:,2)];

top_in = [xs_inbend xs_ext(2:end);ys_inbend ys_inext(2:end)].'*1e-6;
bot_in = [top_in(:,1) -1*top_in(:,2)];

conts = {bot_out;bot_in;top_in;top_out};

end