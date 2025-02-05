clearvars -except wiod_ssm wico_ssm inau_ssm inhu_ssm wcon_ssm ...
    woco_ssm wodu_ssm dgen_ssm
path = '~/Dropbox (MREL)/MATLAB/OO-TechEc/output_data/oossm_out/';
ucs = {'st','lt'};

%preallocate - don't think this buys me much here? ins't correct anyway
% wiod(3,2) = struct();
% wico(3,2) = struct();
% inhu(3,2) = struct();
% inau(3,2) = struct();
% woco(3,2) = struct();
% wcon(3,2) = struct();
% wodu(3,2) = struct();
% dgen(3,2) = struct();

%load data - takes a while
if ~exist('wiod','var') && ~exist('wico','var') && ...
        ~exist('inhu','var') && ~exist('inau','var') && ...
        ~exist('woco','var') && ~exist('wcon','var') && ...
        ~exist('wodu','var') && ~exist('dgen','var')
    tt = tic;
    for l = 1:3
        for uc = 1:2
            disp(['Loading l = ' num2str(l) ' & uc = ' char(ucs(uc)) ...
                ' after ' num2str(toc(tt)/60,2) ' minutes.'])
            wiod_ssm(l,uc) = load([path 'wiod_' ucs{uc} '_' num2str(l)]);
            disp(['      wiod loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            wico_ssm(l,uc) = load([path 'wico_' ucs{uc} '_' num2str(l)]);
            disp(['      wico loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            inhu_ssm(l,uc) = load([path 'inhu_' ucs{uc} '_' num2str(l)]);
            disp(['      inhu loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            inau_ssm(l,uc) = load([path 'inau_' ucs{uc} '_' num2str(l)]);
            disp(['      inau loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            woco_ssm(l,uc) = load([path 'woco_' ucs{uc} '_' num2str(l)]);
            disp(['      woco loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            wcon_ssm(l,uc) = load([path 'wcon_' ucs{uc} '_' num2str(l)]);
            disp(['      wcon loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            wodu_ssm(l,uc) = load([path 'wodu_' ucs{uc} '_' num2str(l)]);
            disp(['      wodu loaded after ' num2str(toc(tt)/60,2) ' mins...'])
            dgen_ssm(l,uc) = load([path 'dgen_' ucs{uc} '_' num2str(l)]);
            disp(['      dgen loaded after ' num2str(toc(tt)/60,2) ' mins...'])
        end
    end
end

%wind
visOoSsmOut('tiv','wiod',wiod_ssm)
visOoSsmOut('tiv','wico',wico_ssm)
visOoSsmOut('tcm','wiod',wiod_ssm)
visOoSsmOut('tcm','wico',wico_ssm)
visOoSsmOut('twf','wiod',wiod_ssm)
visOoSsmOut('twf','wico',wico_ssm)
visOoSsmOut('cis','wiod',wiod_ssm)
visOoSsmOut('cis','wico',wico_ssm)
visOoSsmOut('rsp','wiod',wiod_ssm)
visOoSsmOut('rsp','wico',wico_ssm)
visOoSsmOut('cos','wiod',wiod_ssm) 
visOoSsmOut('cos','wico',wico_ssm)
% visOoSsmOut('tef','wiod',wiod_ssm)
% visOoSsmOut('tef','wico',wico_ssm)
visOoSsmOut('szo','wiod',wiod_ssm)
visOoSsmOut('szo','wico',wico_ssm)
visOoSsmOut('lft','wiod',wiod_ssm) 
visOoSsmOut('lft','wico',wico_ssm)
visOoSsmOut('dtc','wiod',wiod_ssm)
visOoSsmOut('dtc','wico',wico_ssm)
visOoSsmOut('osv','wiod',wiod_ssm)
visOoSsmOut('osv','wico',wico_ssm)
visOoSsmOut('spv','wiod',wiod_ssm)
visOoSsmOut('spv','wico',wico_ssm)
visOoSsmOut('tmt','wiod',wiod_ssm)
visOoSsmOut('tmt','wico',wico_ssm)
visOoSsmOut('eol','wiod',wiod_ssm)
visOoSsmOut('eol','wico',wico_ssm)
visOoSsmOut('dep','wiod',wiod_ssm)
visOoSsmOut('dep','wico',wico_ssm)
visOoSsmOut('bcc','wiod',wiod_ssm)
visOoSsmOut('bcc','wico',wico_ssm)
visOoSsmOut('bhc','wiod',wiod_ssm)
visOoSsmOut('bhc','wico',wico_ssm)
visOoSsmOut('utp','wiod',wiod_ssm)
visOoSsmOut('utp','wico',wico_ssm)
visOoSsmOut('ild','wiod',wiod_ssm)
visOoSsmOut('ild','wico',wico_ssm)
visOoSsmOut('sdr','wiod',wiod_ssm)
visOoSsmOut('sdr','wico',wico_ssm)
%inso
visOoSsmOut('pvd','inau',inau_ssm)
visOoSsmOut('pvd','inhu',inhu_ssm)
visOoSsmOut('pcm','inau',inau_ssm)
visOoSsmOut('pcm','inhu',inhu_ssm)
visOoSsmOut('pwf','inau',inau_ssm)
visOoSsmOut('pwf','inhu',inhu_ssm)
visOoSsmOut('pve','inau',inau_ssm)
visOoSsmOut('pve','inhu',inhu_ssm)
visOoSsmOut('lft','inau',inau_ssm)
visOoSsmOut('lft','inhu',inhu_ssm)
visOoSsmOut('dtc','inau',inau_ssm)
visOoSsmOut('dtc','inhu',inhu_ssm)
visOoSsmOut('osv','inau',inau_ssm)
visOoSsmOut('osv','inhu',inhu_ssm)
visOoSsmOut('spv','inau',inau_ssm)
visOoSsmOut('spv','inhu',inhu_ssm)
visOoSsmOut('tmt','inau',inau_ssm)
visOoSsmOut('tmt','inhu',inhu_ssm)
visOoSsmOut('eol','inau',inau_ssm)
visOoSsmOut('eol','inhu',inhu_ssm)
visOoSsmOut('dep','inau',inau_ssm)
visOoSsmOut('dep','inhu',inhu_ssm)
visOoSsmOut('bcc','inau',inau_ssm)
visOoSsmOut('bcc','inhu',inhu_ssm)
visOoSsmOut('bhc','inau',inau_ssm)
visOoSsmOut('bhc','inhu',inhu_ssm)
visOoSsmOut('utp','inau',inau_ssm)
visOoSsmOut('utp','inhu',inhu_ssm)
visOoSsmOut('ild','inau',inau_ssm)
visOoSsmOut('ild','inhu',inhu_ssm)
visOoSsmOut('sdr','inau',inau_ssm)
visOoSsmOut('sdr','inhu',inhu_ssm)
%wave
visOoSsmOut('wiv','wcon',wcon_ssm)
visOoSsmOut('wiv','woco',woco_ssm)
visOoSsmOut('wiv','wodu',wodu_ssm)
visOoSsmOut('wcm','wcon',wcon_ssm)
visOoSsmOut('wcm','woco',woco_ssm)
visOoSsmOut('wcm','wodu',wodu_ssm)
visOoSsmOut('whl','wcon',wcon_ssm)
visOoSsmOut('whl','woco',woco_ssm)
visOoSsmOut('whl','wodu',wodu_ssm)
visOoSsmOut('ect','wcon',wcon_ssm)
visOoSsmOut('ect','woco',woco_ssm)
visOoSsmOut('ect','wodu',wodu_ssm)
visOoSsmOut('lft','wcon',wcon_ssm)
visOoSsmOut('lft','woco',woco_ssm)
visOoSsmOut('lft','wodu',wodu_ssm)
visOoSsmOut('dtc','wcon',wcon_ssm)
visOoSsmOut('dtc','woco',woco_ssm)
visOoSsmOut('dtc','wodu',wodu_ssm)
visOoSsmOut('osv','wcon',wcon_ssm)
visOoSsmOut('osv','woco',woco_ssm)
visOoSsmOut('osv','wodu',wodu_ssm)
visOoSsmOut('spv','wcon',wcon_ssm)
visOoSsmOut('spv','woco',woco_ssm)
visOoSsmOut('spv','wodu',wodu_ssm)
visOoSsmOut('tmt','wcon',wcon_ssm)
visOoSsmOut('tmt','woco',woco_ssm)
visOoSsmOut('tmt','wodu',wodu_ssm)
visOoSsmOut('eol','wcon',wcon_ssm)
visOoSsmOut('eol','woco',woco_ssm)
visOoSsmOut('eol','wodu',wodu_ssm)
visOoSsmOut('dep','wcon',wcon_ssm)
visOoSsmOut('dep','woco',woco_ssm)
visOoSsmOut('dep','wodu',wodu_ssm)
visOoSsmOut('bcc','wcon',wcon_ssm)
visOoSsmOut('bcc','woco',woco_ssm)
visOoSsmOut('bcc','wodu',wodu_ssm)
visOoSsmOut('bhc','wcon',wcon_ssm)
visOoSsmOut('bhc','woco',woco_ssm)
visOoSsmOut('bhc','wodu',wodu_ssm)
visOoSsmOut('utp','wcon',wcon_ssm)
visOoSsmOut('utp','woco',woco_ssm)
visOoSsmOut('utp','wodu',wodu_ssm)
visOoSsmOut('ild','wcon',wcon_ssm)
visOoSsmOut('ild','woco',woco_ssm)
visOoSsmOut('ild','wodu',wodu_ssm)
visOoSsmOut('sdr','wcon',wcon_ssm)
visOoSsmOut('sdr','woco',woco_ssm)
visOoSsmOut('sdr','wodu',wodu_ssm)
%dgen
visOoSsmOut('giv','dgen',dgen_ssm)
visOoSsmOut('fco','dgen',dgen_ssm)
visOoSsmOut('fca','dgen',dgen_ssm)
visOoSsmOut('fsl','dgen',dgen_ssm)
visOoSsmOut('oci','dgen',dgen_ssm)
visOoSsmOut('gcm','dgen',dgen_ssm)
visOoSsmOut('lft','dgen',dgen_ssm)
visOoSsmOut('dtc','dgen',dgen_ssm)
visOoSsmOut('osv','dgen',dgen_ssm)
visOoSsmOut('spv','dgen',dgen_ssm)
visOoSsmOut('tmt','dgen',dgen_ssm)
visOoSsmOut('eol','dgen',dgen_ssm)
visOoSsmOut('dep','dgen',dgen_ssm)
visOoSsmOut('bcc','dgen',dgen_ssm)
visOoSsmOut('bhc','dgen',dgen_ssm)
visOoSsmOut('utp','dgen',dgen_ssm)
visOoSsmOut('ild','dgen',dgen_ssm)
visOoSsmOut('sdr','dgen',dgen_ssm)





