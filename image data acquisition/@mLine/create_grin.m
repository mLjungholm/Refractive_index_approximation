function grin = create_grin(this, stepNums)
grin = GRIN2d_rotSym(this.refractiveGradient, this.gradientD, stepNums);
end