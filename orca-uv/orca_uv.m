% Plots absorption spectrum from ORCA output files.
% ORCA: https://orcaforum.kofo.mpg.de
% - stick spectra + Gaussian, Lorentzian & Pseudo-Voigt line broadening
% - nm or cm^-1^ scale
% - several options
% 
% 1. Check the options.
% 2. Run the script with F5. 
% 3. Open an ORCA output file.
% 4. Spectrum will be saved as 'abs-spectrum.png' in the same folder
%    - You need write permissions!
%    - A previous 'abs-spectrum.png' file will be overwritten!
%

clear all;

% options
tedp = 1;                                % 1 for TRANSITION ELECTRIC DIPOLE MOMENTS, 0 for TRANSITION VELOCITY DIPOLE MOMENTS
w_nm = 20;                               % peak broadening for nm scale
w_wn = 2000;                             % peak broadening for wavenumber scale
start_x_nm = 0;                          % start spectrum at ... nm
start_x_wn = 0;                          % start spectrum at ... cm^-1^ 
add_to_max_x_nm = 0;                     % expand x plus ... nm (-... decrease)
add_to_max_x_wn = 12000;                 % expand x plus ... cm^-1^ (-... decrease)
nm = 1;                                  % 1 for wavelength (nm), 0 for wavenumber (cm^-1^)
gaussian_ls = 1;                         % Gaussian line shape
lorentzian_ls = 0;                       % Lorentzian line shape
pvoigt_ls = 0;                           % Pseudo-Voigt line shape
gauss_area = 1;                          % area plot for Gaussian - else line
lorentz_area = 0;                        % area plot for Lorentzian - else line
pvoigt_area = 0;                         % area plot for Pseudo-Voigt - else line
gauss_single = 0;                        % plot Gaussian for every single peak
lorentz_single = 0;                      % plot Lorentzian for every single peak
pvoigt_single = 0;                       % plot Pseudo-Voigt for every single peak
hiwn_to_lown = 1;                        % spectrum starts from high wavenumber
peak_detection = 1;                      % turn peak detection on
peak_treshold = 0.001;                   % peak detection treshold
peak_font_size = 7;                      % font size for detected peaks
area_color = [0.3 0.3 0.3];              % area color
stick_color = [0 0 0];                   % stick color
spectrum_title = 'Absorption spectrum';  % spectrum title
resolution = 300;                        % resolution of the picture in dpi 
% end options

extracted_uv_data={};

% open file and get ID
[orca_out_file, path] = uigetfile('.out','Select an ORCA output file');
if isequal(orca_out_file,0)
   error('Select an ORCA output file!');
   return;
else
   orca_out_file_ID=fopen(fullfile(path,orca_out_file),'r');
end

% look for string 'ABSORPTION SPECTRUM' in orca.out
while ~feof(orca_out_file_ID)
    line=fgetl(orca_out_file_ID);
    if tedp
        is_uv_spectrum=strfind(line,'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS');
    else
        is_uv_spectrum=strfind(line,'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS');
    end
    if is_uv_spectrum
        break
    end
end

if isempty(is_uv_spectrum)
    error('Absorption spectrum not found!') % error mesage > end prgm
else
    % extract lines with IR data
    % skip next 4 lines
    for nlines_2_skip = 1:4
        fgetl(orca_out_file_ID);
    end
    % read IR data 
    line = fgetl(orca_out_file_ID);
    while ~isempty(line)
       extracted_uv_data_line=textscan(line,'%s %f %f %f %*[^\n]');
       extracted_uv_data_temp=[extracted_uv_data_line{:}];
       extracted_uv_data(end+1,:) = extracted_uv_data_temp;
       line = fgetl(orca_out_file_ID);
   end
end
%
fclose(orca_out_file_ID);

x_wn = [extracted_uv_data{:,2}]; % x (wavenumbers) from orca.out
x_nm = [extracted_uv_data{:,3}]; % x (nm) from orca.out
y = [extracted_uv_data{:,4}];    % y (intesities) from orca.out

if nm
    x = x_nm;
    w = w_nm;
    add_to_max_x = add_to_max_x_nm;
    start_x = start_x_nm;
    box_start_x = 0;
else
    x = x_wn;
    w = w_wn;
    add_to_max_x = add_to_max_x_wn;
    start_x = start_x_wn;
    if hiwn_to_lown
        box_start_x = start_x_wn;
    else
        box_start_x = 0;
    end
end

% x_0 should not be lt 0 or exceed x_max 
if start_x >= max(x) || start_x == 0
    start_x = 0;    
end

figure('Name','Absorption spectrum');

hold on; % several plots in one figure

% high wavenumbers to low wavenumbers
if hiwn_to_lown && nm ~= 1
    set(gca,'XDir','reverse'); 
end

xlim([start_x max(x)+add_to_max_x]); % x plot limits
%z=0:1:max(x)+add_to_max_x;  % x limits for gaussian, lorentzian or pseudo-voigt
z=1:max(x)+add_to_max_x;     % x limits for gaussian, lorentzian or pseudo-voigt

% plot gaussian for every single peak
if gauss_single
    for i = 1 : numel(x)
        ar=area(gauss(y(i),z,x(i),w));
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    end
end

% plot lorentzian for every single peak
if lorentz_single
    for i = 1 : numel(x)
        ar=area(lorentz(y(i),z,x(i),w));
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    end
end

% plot pseudo voigt for every single peak
if pvoigt_single
    for i = 1 : numel(x)
        ar=area(pvoigt(y(i),z,x(i),w));
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    end
end

% plot gaussian
if gaussian_ls
    gauss_sum=sum(gauss(y(:),z,x(:),w));
    yhs_mg = max(gauss_sum)*0.15;  % headspace in y (to avoid labels crossed by lines)
    ylim([0 max(gauss_sum)+yhs_mg]);
    % area or line plot
    if gauss_area
        ar=area(gauss_sum);
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    else
        pl=plot(gauss_sum);
        pl.Color=[0 0 0];
    end
    % begin peak detection
    if peak_detection
        label_dist=max(gauss_sum)*0.01; % peak labeling from top
        [maxtab, mintab] = peakdet(gauss_sum,0.0001,z);
        for i = 1 : numel(maxtab(:,2))
            if maxtab(i,2) > peak_treshold % treshold
                tx=text(maxtab(i,1), maxtab(i,2)+label_dist,num2str(round(maxtab(i,1),0)));
                set(tx,'Rotation',90); % for absorption
                set(tx,'FontSize',peak_font_size);
            end
        end
    end
    % end peak detection
end

% plot lorentzian
if lorentzian_ls
    lorentz_sum=sum(lorentz(y(:),z,x(:),w));
    yhs_ml = max(lorentz_sum)*0.15;  % headspace in y (to avoid labels crossed by lines)
    ylim([0 max(lorentz_sum)+yhs_ml]);
    % area or line plot
    if lorentz_area
        ar=area(lorentz_sum);
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    else
        pl=plot(lorentz_sum);
        pl.Color=[0 0 0];
    end
    % begin peak detection
    if peak_detection
        label_dist=max(lorentz_sum)*0.01; % peak labeling from top
        [maxtab, mintab] = peakdet(lorentz_sum,0.0001,z);
        for i = 1 : numel(maxtab(:,2))
            if maxtab(i,2) > peak_treshold % treshold
                tx=text(maxtab(i,1), maxtab(i,2)+label_dist,num2str(round(maxtab(i,1),0)));
                set(tx,'Rotation',90); % for absorption
                set(tx,'FontSize',peak_font_size);
            end
        end
    end
    % end peak detection
end

% plot pseudo-voigt
if pvoigt_ls
    pvoigt_sum=sum(pvoigt(y(:),z,x(:),w));
    yhs_pv = max(pvoigt_sum)*0.15;  % headspace in y (to avoid labels crossed by lines)
    ylim([0 max(pvoigt_sum)+yhs_pv]);
    % area or line plot
    if pvoigt_area
        ar=area(pvoigt_sum);
        ar.FaceColor = area_color;
        ar.EdgeColor = [0 0 0];
        ar.FaceAlpha = 0.2;
        ar.EdgeAlpha = 0.6;
    else
        pl=plot(pvoigt_sum);
        pl.Color=[0 0 0];
    end
    % begin peak detection
    if peak_detection
        label_dist=max(pvoigt_sum)*0.01; % peak labeling from top
        [maxtab, mintab] = peakdet(pvoigt_sum,0.0001,z);
        for i = 1 : numel(maxtab(:,2))
            if maxtab(i,2) > peak_treshold % treshold
                tx=text(maxtab(i,1), maxtab(i,2)+label_dist,num2str(round(maxtab(i,1),0)));
                set(tx,'Rotation',90); % for absorption
                set(tx,'FontSize',peak_font_size);
            end
        end
    end
    % end peak detection
end

% plot stick spectrum
st=stem(x,y);
st.Marker='none';
st.Color=stick_color;

% plot a box 
if gaussian_ls && lorentzian_ls == 0 && pvoigt_ls == 0
    rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(gauss_sum)+yhs_mg]);
    ylim([0 max(gauss_sum)+yhs_mg]);
end
if lorentzian_ls && gaussian_ls == 0 && pvoigt_ls == 0
    rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(lorentz_sum)+yhs_ml]);
    ylim([0 max(lorentz_sum)+yhs_ml]);
end
if pvoigt_ls && gaussian_ls == 0 && lorentzian_ls == 0
    rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(pvoigt_sum)+yhs_pv]);
    ylim([0 max(pvoigt_sum)+yhs_pv]);
end
if lorentzian_ls && gaussian_ls && pvoigt_ls == 0
    if max(gauss_sum) > max(lorentz_sum)
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(gauss_sum)+yhs_mg]);
        ylim([0 max(gauss_sum)+yhs_mg]);
    else
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(lorentz_sum)+yhs_ml]);
        ylim([0 max(lorentz_sum)+yhs_ml]);
    end
end
if lorentzian_ls && pvoigt_ls && gaussian_ls == 0
    if max(pvoigt_sum) > max(lorentz_sum)
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(pvoigt_sum)+yhs_pv]);
        ylim([0 max(pvoigt_sum)+yhs_pv]);
    else
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(lorentz_sum)+yhs_ml]);
        ylim([0 max(lorentz_sum)+yhs_ml]);
    end
end
if pvoigt_ls && gaussian_ls && lorentzian_ls == 0
    if max(gauss_sum) > max(pvoigt_sum)
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(gauss_sum)+yhs_mg]);
         ylim([0 max(gauss_sum)+yhs_mg]);
    else
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(pvoigt_sum)+yhs_pv]);
        ylim([0 max(pvoigt_sum)+yhs_pv]);
    end
end
if pvoigt_ls && gaussian_ls && lorentzian_ls
    if max(gauss_sum) > max(pvoigt_sum) && max(gauss_sum) > max(lorentz_sum)
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(gauss_sum)+yhs_mg]);
        ylim([0 max(gauss_sum)+yhs_mg]);
    elseif max(pvoigt_sum) > max(gauss_sum) && max(pvoigt_sum) > max(lorentz_sum)
        rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(pvoigt_sum)+yhs_pv]);
        ylim([0 max(pvoigt_sum)+yhs_pv]);
    elseif max(lorentz_sum) > max(gauss_sum) && max(lorentz_sum) > max(pvoigt_sum)
        re=rectangle('Position',[box_start_x 0 max(x)+add_to_max_x max(lorentz_sum)+yhs_ml]);
        ylim([0 max(lorentz_sum)+yhs_ml]);
    end
end

% some plot options
box off;
%grid on;
%grid minor;
set(gca,'XMinorTick','on');
title(spectrum_title);
if nm
    xlabel(['wavelength / nm']);
else
    xlabel(['wavenumber / cm^{' char(8211) '1}']);
end
ax = gca;
ax.XAxis.Exponent = 0;
ylabel('intensity');
set(gca,'ytick',[])
set(gca,'TickDir','out');
%set(gca,'linewidth',1);

hold off;

% save gfx as png; resolution in dpi
% same folder as orca.out
print([path 'abs-spectrum'],'-dpng',['-r' num2str(resolution)]);

% calculation of the Gaussian line shape
% a = amplitude (max y, intensity)
% x = position
% m = maximum/meadian (stick position in x, wavelength or wavenumber)
% w = line width, FWHM
function n = gauss(a,x,m,w)
    % calculation of Gaussian
    %n=a.*1/(s*sqrt(2*pi)).*exp(-((x-m).^2/(2*s.^2)));
    n=a.*exp(-(log(2).*((m-x)/w).^2));
end

% calculation of the Lorentzian line shape
% a = amplitude (max y, intensity)
% x = position
% m = maximum/meadian (stick position in x, wavelength or wavenumber)
% w = line width, FWHM
function n = lorentz(a,x,m,w)
   % calculation of Lorentzian
   n=a.*1./(1+((m-x)/w).^2);
end

% calculation of the Pseudo-Voigt line shape
% a = amplitude (max y, intensity)
% x = position
% m = maximum/meadian (stick position in x, wavelength or wavenumber)
% w = line width, FWHM
function n = pvoigt(a,x,m,w)
  % calculation of Pseudo-Voigt: pvoigt = h*lorentz + (1-h)*gauss; 0 < h < 1
  h=0.5;
  n=h.*(a.*1./(1+((m-x)/w).^2))+(1-h).*(a.*exp(-(log(2).*((m-x)/w).^2)));
end

function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end
end

% BSD 3-Clause License
% 
% Copyright (c) 2019, Sebastian Dechert
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
