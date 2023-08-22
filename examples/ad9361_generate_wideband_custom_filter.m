% Requires the Communications Toolbox Support Package for Xilinx Zynq-Based
% Radio

% Please read https://github.com/analogdevicesinc/ad936x-filter-wizard/blob/master/cook_input.m
% for more information on the filter design parameters

% The development of this code is based on the analysis of design_filter.m and
% more specifically internal_design_filter.m
% https://github.com/analogdevicesinc/ad936x-filter-wizard/blob/master/design_filter.m
% https://github.com/analogdevicesinc/ad936x-filter-wizard/blob/master/internal_design_filter.m


% Create a partial structure for the filter design: 
fdp.Rdata   = 61.44e6;  % Output data rate in Hz
fdp.Fpass   = 56e6;   % FIR passband edge in Hz
fdp.Fstop   = 58e6;   % FIR stopband edge in Hz
fdp.Fcutoff = 56e6;     % Cutoff frequency in Hz
fdp.Astop   = 30;       % Stopband attenuation in dB

% Complete the filter design structure
fdp = cook_input(fdp)

% Perform the filter design
fir_struct = design_filter(fdp)

% Plot the filter response
fvtool(fir_struct.firtaps, 'Fs', fdp.Rdata)


export_filter_flag = false;

if export_filter_flag
    
% Export the filter coefficients to a file
% based on save2coefficients_Callback in AD9361_Filter_Wizard.m

fid = fopen('61_44_MHz_custom_fir.ftr','w');

fprintf(fid, '# Generated with AD9361 Filter Design Wizard %s\r\n', get_version);
fprintf(fid, '# MATLAB %s, %s\r\n', version(), datestr(now()));
fprintf(fid, '# Inputs:\r\n');

PLL_rate = fir_struct.PLL_rate;

rx_FIR_rate = fir_struct.Rdata * fir_struct.FIR;
rx_HB1_rate = rx_FIR_rate * fir_struct.HB1;
rx_HB2_rate = rx_HB1_rate * fir_struct.HB2;
rx_HB3_rate = rx_HB2_rate * fir_struct.HB3;

fprintf(fid, '# Data Sample Frequency = %.0f Hz\r\n', fir_struct.Rdata);
fprintf(fid, '# RX Phase equalization = %f ns\r\n', fir_struct.phEQ);
fprintf(fid, 'RX 3 GAIN %d DEC %d\r\n', fir_struct.gain, fir_struct.FIR);
fprintf(fid, 'RRX %.0f %.0f %.0f %.0f %.0f %.0f\r\n', PLL_rate, rx_HB3_rate, rx_HB2_rate, rx_HB1_rate, rx_FIR_rate, fir_struct.Rdata);
fprintf(fid, 'BWRX %.0f\r\n', fir_struct.RFbw);

% concat and transform Rx and Tx coefficient matrices for output
coefficients = flip(rot90(int16(fir_struct.firtaps)));

% output all non-zero coefficients since they're padded to 128 with zeros
for i = 1:fir_struct.nfirtaps
    fprintf(fid, '%d\r\n', coefficients(i));
end

fclose(fid);

end