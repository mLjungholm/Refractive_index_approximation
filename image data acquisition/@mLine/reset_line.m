function reset_line(this,options)
arguments
    this
    options {mustBeMember(options,{'all','peaks','edges','phase'})} = 'all'
end
switch options
    case 'peaks'
        this.gaussPks = {};
        this.sgolayPks = {};
        this.pks = {};
        this.PD = [];
        this.PL = [];
        this.PS = [];
        this.PP = [];
        this.PV = [];
        this.L = {};
        this.S = {};
        this.reset_line('phase')
        
    case 'phase'
        this.phasePoints = [];
        
    case 'edges'
        this.leftEdge = nan;
        this.rigthEdge = nan;
        this.leftPhaseMax = nan;
        this.rightPhaseMin = nan;
end





end