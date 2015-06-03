% im = colorborder(im, bord, color)
%
% This function puts a border of width 'bord' around the image 'im'
% using the color [r g b].
%
% Example:
% >> im = colorborder(im, 5, [1 0 0]);  % puts a red border
function im = colorborder(im, bord, color)

im = padarray(im, [bord bord 0], 0);

for chan=1:3,
  im(1:bord,         :, chan) = color(chan);
  im(end-bord+1:end, :, chan) = color(chan);
  im(:, 1:bord,         chan) = color(chan);
  im(:, end-bord+1:end, chan) = color(chan);
end