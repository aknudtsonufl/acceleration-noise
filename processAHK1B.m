function Temps = processAHK1B(fname)
    % Load temperature data only from HK1B files

    % open file and load data
    fid = fopen(fname);
    disp(fid)
    dd = textscan(fid,'%f%f%c%c%s%s%f%f%f%f%f%s\n','headerlines', 268);
    fclose(fid);

    % launch time
    tLaunch = time('2018-05-22 17:37:19.000');

    % bit mapping
    % bit/ind 00/32/24 = TFEEU_IF
    % bit/ind 01/31/23 = TFEEU_REF
    % bit/ind 02/30/22 = TFEEU_X
    % bit/ind 03/29/21 = TFEEU_YZ
    % bit/ind 04/28/20 = analog_GND
    % bit/ind 05/27/19 = +3.3V
    % bit/ind 06/26/18 = Vp
    % bit/ind 07/25/-- = MES_Vd
    % bit/ind 08/24/17 = MES_DET_X1
    % bit/ind 09/23/16 = MES_DET_X2
    % bit/ind 10/22/15 = MES_DET_X3
    % bit/ind 11/21/14 = MES_DET_Y1
    % bit/ind 12/20/13 = MES_DET_Y2
    % bit/ind 13/19/12 = MES_DET_Z1
    % bit/ind 14/18/11 = TSU_Y+
    % bit/ind 15/17/10 = TICUN
    % bit/ind 16/16/09 = TSU_Y-
    % bit/ind 17/15/08 = TSU_Z+
    % bit/ind 18/14/07 = TSU_Z-
    % bit/ind 19/13/06 = +5V
    % bit/ind 20/12/05 = TICUR
    % bit/ind 21/11/04 = +15V
    % bit/ind 22/10/03 = -15V
    % bit/ind 23/09/02 = +48V
    % bit/ind 24/08/01 = -48V

    AHKLabels = {'TFEEU_IF', 'TFEEU_REF', 'TFEEU_X', 'TFEEU_YZ', 'analog_GND', ...
      '+3.3V', 'Vp', '', 'MES_DET_X1', 'MES_DET_X2', 'MES_DET_X3', 'MES_DET_Y1', ...
      'MES_DET_Y2', 'MES_DET_Z1', 'TSU_Y+', 'TICUN', 'TSU_Y-', 'TSU_Z+', 'TSU_Z-', ...
      '+5V', 'TICUR', '+15V', '-15V', '+48V', '-48V'};
    AHKUnits = {'degC', 'degC', 'degC', 'degC', 'V', ...
      'V', 'V', '','m', 'm', 'm', 'm', ...
      'm', 'm', 'degC', 'degC', 'degC', 'degC', 'degC', ...
      'V', 'degC', 'V', 'V', 'V', 'V'};

    % spacecraft reference
    disp(dd)
    sc = dd{4}(1);

    % data container
    HK = [];

    % data fields
    d1 = dd{7};
    d2 = dd{8};

    % time in seconds
    ts = dd{1};
    % factional seconds
    tf = dd{2}/1e6;
    % total time
    tt = ts+tf;

    count = 1;


    indName = [12 14:18 29:32];
    % Loop through desired HK values in bit-map
    for jj = [12 14:18 29:32]%[8:24 26:32]
      clear ind
      % Loop through all data, find the indices where the desired bit set to 1
      for ii = 1:length(dd{6})
        ind(ii) = strcmp(dd{6}{ii}(jj),'1');
      end
      % TM polarization is always present. If desired HK is before it in bit
      % list, desired HK is 2nd value, otherwise the 1st
      if jj<25
        HK = [HK ao(plist('type', 'tsdata', 'xvals', tt(ind), 'yvals', d2(ind), ...
          'yunits', AHKUnits{length(AHKUnits)+1-count}, 't0', '2000-Jan-01 00:00:00'))];
    %     HK(count).setName([AHKLabels{length(AHKLabels)+8-jj} sprintf('_%c',sc)]);
      else
        HK = [HK ao(plist('type', 'tsdata', 'xvals', tt(ind), 'yvals', d1(ind), ...
          'yunits', AHKUnits{length(AHKUnits)+1-count}, 't0', '2000-Jan-01 00:00:00'))];
    %     HK(count).setName([AHKLabels{length(AHKLabels)+8-jj} sprintf('_%c',sc)]);
      end
      count = count + 1;
    end

    HK.setReferenceTime(tLaunch);

    % Fix the sampling frequency, interpolate to same grip
    Temps = consolidate(HK, plist('fs', 0.2));

    for ii = 1:numel(indName)
      Temps(ii).setName([AHKLabels{length(AHKLabels)+8-indName(ii)} sprintf('_%c',sc)]);
    end
    
end