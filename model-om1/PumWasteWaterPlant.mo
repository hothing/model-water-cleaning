model PumWasteWaterPlant

  type MassConcentration = Real (
    final quantity="MassConcentration",
    final unit="mg/l",
    final min=0.0);

  type Alkalinity = Real (
    final quantity="Alkalinity",
    final unit="mmol/l",
    final min=0.0);

  type VolumeFlowRate = Real (final quantity="VolumeFlowRate", final unit= "m3/d");

  connector FlowIn
    flow VolumeFlowRate Q;
    MassConcentration Si;
    MassConcentration Ss;
    MassConcentration Xi;
    MassConcentration Xs;
    MassConcentration Xbh;
    MassConcentration Xba;
    MassConcentration Xp;
    MassConcentration So;
    MassConcentration Sno;
    MassConcentration Snh;
    MassConcentration Snd;
    MassConcentration Xnd;
    Alkalinity Salk;
  annotation(
      Icon(graphics = {Polygon(fillColor = {52, 101, 164}, fillPattern = FillPattern.Solid, points = {{-60, 60}, {60, 60}, {60, 0}, {60, -60}, {-60, -60}, {-20, 0}, {-60, 60}, {-20, 0}, {-60, 60}})}));
  end FlowIn;
  
  connector FlowOut
    flow VolumeFlowRate Q;
    MassConcentration Si;
    MassConcentration Ss;
    MassConcentration Xi;
    MassConcentration Xs;
    MassConcentration Xbh;
    MassConcentration Xba;
    MassConcentration Xp;
    MassConcentration So;
    MassConcentration Sno;
    MassConcentration Snh;
    MassConcentration Snd;
    MassConcentration Xnd;
    Alkalinity Salk;
  annotation(
      Icon(graphics = {Polygon(origin = {-10, 0}, fillColor = {115, 210, 22}, fillPattern = FillPattern.Solid, points = {{-50, 60}, {30, 60}, {70, 0}, {30, -60}, {-50, -60}, {-50, 0}, {-50, 60}})}));
  end FlowOut;
  
  model ASM1
    //inputs
    PumWasteWaterPlant.FlowIn flowIn annotation(
      Placement(visible = true, transformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    PumWasteWaterPlant.FlowOut flowOut annotation(
      Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput T annotation(
      Placement(visible = true, transformation(origin = {-108, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput aeration annotation(
      Placement(visible = true, transformation(origin = {-108, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
    // parameters
    // tank specific parameters
    parameter Modelica.SIunits.Volume V=1000 "Volume of denitrification tank";
    
    parameter Real mu_h_T = 4.0 "Maximum heterotrophic growth rate at T=15 deg C [day^-1]";
    parameter Real b_h_T = 0.28 "Heterotrophic decay rate at T=15 deg C [day^-1]";
    parameter Real mu_a_T = 0.5 "Maximum autotrophic growth rate at T=15 deg C[day^-1]";
    parameter Real b_a_T = 0.1 "Autotrophic decay rate at T=15 deg C [day^-1]";
    parameter Real k_a_T = 0.06 "Ammonification rate at T=15 deg C [m3/(g COD day)]";
    parameter Real k_h_T = 1.75 "Maximum specific hydrolysis rate at T=15 deg C [g Xs/(g Xbh COD day)]";
    parameter Real K_x_T = 0.0175 "Half-saturation (hydrolysis) at T=15 deg C [g Xs/(g Xbh COD)]";
    parameter Real K_nh = 1.0 "Half-saturation (auto. growth) [g NH-N/m3]";
    parameter Real K_s = 20.0 "Half-saturation (hetero. growth) [g COD/m3]";
    parameter Real K_oh = 0.2 "Half-saturation (hetero. oxygen) [g O/m3]";
    parameter Real K_no = 0.5 "Half-saturation (nitrate) [g NO-N/m3]";
    parameter Real K_oa = 0.4 "Half-saturation (auto. oxygen) [g O/m3]";
    parameter Real ny_g = 0.8 "Anoxic growth rate correction factor [-]";
    parameter Real ny_h = 0.4 "Anoxic hydrolysis rate correction factor [-]";
    // Stoichiometric parameters based on the original ASM1 publication//
    parameter Real Y_h=0.67
      "Heterotrophic Yield [g Xbh COD formed/(g COD utilised)]";
    parameter Real Y_a=0.24
      "Autotrophic Yield [g Xba COD formed/(g N utilised)]";
    parameter Real f_p=0.08 "Fraction of biomass to particulate products [-]";
    parameter Real i_xb=0.086 "Fraction nitrogen in biomass [g N/(g COD)]";
    parameter Real i_xp=0.06
      "Fraction nitrogen in particulate products [g N/(g COD)]";
    
    // parameters based on the original ASM1 publication based on 15 deg C
    Real mu_h "Maximum heterotrophic growth rate f(T) [day^-1]";
    Real b_h "Heterotrophic decay rate f(T) [day^-1]";
    Real mu_a "Maximum autotrophic growth rate f(T) [day^-1]";
    //Real K_nh "Half-saturation (auto. growth) f(T) [g NH-N/m3]";
    Real b_a "Autotrophic decay rate f(T) [day^-1]";
    Real k_a "Ammonification rate f(T) [m3/(g COD day)]";
    Real k_h "Maximum specific hydrolysis rate f(T) [g Xs/(g Xbh COD day)]";
    Real K_x "Half-saturation (hydrolysis) f(T) [g Xs/(g Xbh COD)]";
    
    // internal variables
    Real p1;
    Real p2;
    Real p3;
    Real p4;
    Real p5;
    Real p6;
    Real p7;
    Real p8;
    Real r1;
    Real r2;
    Real r3;
    Real r4;
    Real r5;
    Real r6;
    Real r7;
    Real r8;
    Real r9;
    Real r10;
    Real r11;
    Real r12;
    Real r13;
  
    // process variables (state variables)
    MassConcentration Si(fixed=true) "Soluble inert organic matter";
    MassConcentration Ss(fixed=true) "Readily biodegradable substrate";
    MassConcentration Xi(fixed=true) "Particulate inert organic matter";
    MassConcentration Xs(fixed=true) "Slowly biodegradable substrate";
    MassConcentration Xbh(fixed=true) "Active heterotrophic biomass";
    MassConcentration Xba(fixed=true) "Active autotrophic biomass";
    MassConcentration Xp(fixed=true) "Particulate products from biomass decay";
    MassConcentration So(fixed=true) "Dissolved oxygen";
    MassConcentration Sno(fixed=true) "Nitrate and nitrite nitrogen";
    MassConcentration Snh(fixed=true) "Ammonium nitrogen";
    MassConcentration Snd(fixed=true) "Soluble biodegradable organic nitrogen";
    MassConcentration Xnd(fixed=true) "Particulate biodegradable organic nitrogen";
    Alkalinity Salk(fixed=true) "Alkalinity";
  
  equation

    // Temperature dependent Kinetic parameters based on 15 deg C //
  // may be adapted to 10 or 20 deg C
    mu_h = mu_h_T * exp(0.069 * (T - 15));
    b_h = b_h_T * exp(0.069 * (T - 15));
    mu_a = mu_a_T * exp(0.098 * (T - 15));
    b_a = b_a_T * exp(0.08 * (T - 15));
    k_a = k_a_T * exp(0.069 * (T - 15));
    k_h = k_h_T * exp(0.11 * (T - 15));
    K_x = K_x_T * exp(0.11 * (T - 15));
    
    // Process Rates
    p1 = mu_h * (Ss / (K_s + Ss)) * (So / (K_oh + So)) * Xbh;
    p2 = mu_h * (Ss / (K_s + Ss)) * (K_oh / (K_oh + So)) * (Sno / (K_no + Sno)) * ny_g * Xbh;
    p3 = mu_a * (Snh / (K_nh + Snh)) * (So / (K_oa + So)) * Xba;
    p4 = b_h * Xbh;
    p5 = b_a * Xba;
    p6 = k_a * Snd * Xbh;
    p7 = k_h * (Xs / Xbh / (K_x + Xs / Xbh)) * (So / (K_oh + So) + ny_h * (K_oh / (K_oh + So)) * (Sno / (K_no + Sno))) * Xbh;
    p8 = p7 * Xnd / Xs;

    // biochemical reactions
    r1 = 0;
    r2 = ((-p1) - p2) / Y_h + p7;
    r3 = 0;
    r4 = (1 - f_p) * (p4 + p5) - p7;
    r5 = p1 + p2 - p4;
    r6 = p3 - p5;
    r7 = f_p * (p4 + p5);
    r8 = (-(1 - Y_h) / Y_h * p1) - (4.57 - Y_a) / Y_a * p3;
    r9 = (-(1 - Y_h) / (2.86 * Y_h) * p2) + p3 / Y_a;
    r10 = (-i_xb * (p1 + p2)) - (i_xb + 1 / Y_a) * p3 + p6;
    r11 = (-p6) + p8;
    r12 = (i_xb - f_p * i_xp) * (p4 + p5) - p8;
    r13 = (-i_xb / 14 * p1) + ((1 - Y_h) / (14 * 2.86 * Y_h) - i_xb / 14) * p2 - (i_xb / 14 + 1 / (7 * Y_a)) * p3 + p6 / 14;
    
    // ODE systems
    der(Si) = flowIn.Si + r1;
    der(Ss) = flowIn.Ss + r2;
    der(Xi) = flowIn.Xi + r3;
    der(Xs) = flowIn.Xs + r4;
    der(Xbh) = flowIn.Xbh + r5;
    der(Xba) = flowIn.Xba + r6;
    der(Xp) = flowIn.Xp + r7;
    der(So) = flowIn.So + r8 + aeration;
    der(Sno) = flowIn.Sno + r9;
    der(Snh) = flowIn.Snh + r10;
    der(Snd) = flowIn.Snd + r11;
    der(Xnd) = flowIn.Xnd + r12;
    der(Salk) = flowIn.Salk + r13;
    
    // Outputs
  
    flowOut.Q + flowIn.Q = 0;
    
    flowOut.Si = Si;
    flowOut.Ss = Ss;
    flowOut.Xi = Xi;
    flowOut.Xs = Xs;
    flowOut.Xbh = Xbh;
    flowOut.Xba = Xba;
    flowOut.Xp = Xp;
    flowOut.So = So;
    flowOut.Sno = Sno;
    flowOut.Snh = Snh;
    flowOut.Snd = Snd;
    flowOut.Xnd = Xnd;
    flowOut.Salk = Salk;
    
  annotation(
      Icon(graphics = {Rectangle(origin = {0, 1}, fillColor = {204, 0, 0}, fillPattern = FillPattern.Backward, lineThickness = 1, extent = {{-100, 99}, {100, -99}})}));end ASM1;

  model Test
    PumWasteWaterPlant.ASM1 asm11 annotation(
      Placement(visible = true, transformation(origin = {30, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    PumWasteWaterPlant.WaterSource waterSource1(Salk = 1, Si = 1, Snd = 1, Snh = 1, Sno = 1, So = 1, Ss = 1, Xba = 1, Xbh = 1, Xi = 1, Xnd = 1, Xp = 1, Xs = 1, gFlow = 1)  annotation(
      Placement(visible = true, transformation(origin = {-52, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant Tamb(k = 30)  annotation(
      Placement(visible = true, transformation(origin = {-24, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant flowAir(k = 1)  annotation(
      Placement(visible = true, transformation(origin = { -24, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step stepInflow annotation(
      Placement(visible = true, transformation(origin = {-82, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  PumWasteWaterPlant.WaterSink waterSink1 annotation(
      Placement(visible = true, transformation(origin = {70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
  connect(asm11.flowOut, waterSink1.flowIn) annotation(
      Line(points = {{41, 16}, {36, 16}, {36, 17}, {63, 17}}));
  connect(waterSource1.flowOut, asm11.flowIn) annotation(
      Line(points = {{-43, 15}, {19, 15}, {19, 16}}));
  connect(Tamb.y, asm11.T) annotation(
      Line(points = {{-12, -2}, {-8, -2}, {-8, 10}, {20, 10}}, color = {0, 0, 127}));
  connect(flowAir.y, asm11.aeration) annotation(
      Line(points = {{-12, -36}, {-4, -36}, {-4, 6}, {20, 6}}, color = {0, 0, 127}));
    connect(stepInflow.y, waterSource1.u) annotation(
      Line(points = {{-70, 24}, {-66, 24}, {-66, 22}, {-58, 22}, {-58, 22}}, color = {0, 0, 127}));
  end Test;


  model WaterSource
  
    parameter Real gFlow;
    parameter MassConcentration Si;
    parameter MassConcentration Ss;
    parameter MassConcentration Xi;
    parameter MassConcentration Xs;
    parameter MassConcentration Xbh;
    parameter MassConcentration Xba;
    parameter MassConcentration Xp;
    parameter MassConcentration So;
    parameter MassConcentration Sno;
    parameter MassConcentration Snh;
    parameter MassConcentration Snd;
    parameter MassConcentration Xnd;
    parameter Alkalinity Salk;
    
    PumWasteWaterPlant.FlowOut flowOut annotation(
      Placement(visible = true, transformation(origin = {94, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {94, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput u annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-52, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
  
    flowOut.Q = -gFlow * u;
    flowOut.Si = Si;
    flowOut.Ss = Ss;
    flowOut.Xi = Xi;
    flowOut.Xs = Xs;
    flowOut.Xbh = Xbh;
    flowOut.Xba = Xba;
    flowOut.Xp = Xp;
    flowOut.So = So;
    flowOut.Sno = Sno;
    flowOut.Snh = Snh;
    flowOut.Snd = Snd;
    flowOut.Xnd = Xnd;
    flowOut.Salk = Salk;
    
    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={
          Ellipse(
            extent={{-54,54},{56,-54}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-8,-54},{-14,-52},{-24,-48},{-32,-44},{-36,-40},{-42,-34},{
                -48,-26},{-50,-20},{52,-20},{50,-26},{46,-32},{42,-36},{38,-40},{
                34,-44},{30,-46},{26,-48},{22,-50},{16,-52},{10,-54},{4,-54},{0,
                -54},{-8,-54}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-54,54},{56,-54}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Rectangle(
            extent={{-4,-52},{4,-74}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{4,-74},{88,-68}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-4,-54},{-4,-74},{88,-74}},
            thickness=0.5),
          Line(
            points={{4,-54},{4,-68},{88,-68}},
            thickness=0.5)}),
      Documentation(info="This component is used to feed an ASM1 wwtp model with flow data from measurement
  when e.g. concentration is measured after the primary clarifier.
  
  The dimension of InPort is 1.
  
    1 volumeflowrate Q of incoming wastewater [m3/d]"),
  Icon(graphics={
          Ellipse(
            extent={{-54,54},{56,-54}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-8,-54},{-14,-52},{-24,-48},{-32,-44},{-36,-40},{-42,-34},{
                -48,-26},{-50,-20},{52,-20},{50,-26},{46,-32},{42,-36},{38,-40},{
                34,-44},{30,-46},{26,-48},{22,-50},{16,-52},{10,-54},{4,-54},{0,
                -54},{-8,-54}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-54,54},{56,-54}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Rectangle(
            extent={{-4,-52},{4,-74}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{4,-74},{88,-68}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,95,191},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-4,-54},{-4,-74},{88,-74}},
            thickness=0.5),
          Line(
            points={{4,-54},{4,-68},{88,-68}},
            thickness=0.5)}));
  
  end WaterSource;

  model WaterSink
    PumWasteWaterPlant.FlowIn flowIn annotation(
      Placement(visible = true, transformation(origin = {-76, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-74, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real qi;
  equation
    qi = -flowIn.Q;
    annotation(
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Ellipse(extent = {{-54, 54}, {56, -54}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{-8, -54}, {-14, -52}, {-24, -48}, {-32, -44}, {-36, -40}, {-42, -34}, {-48, -26}, {-50, -20}, {52, -20}, {50, -26}, {46, -32}, {42, -36}, {38, -40}, {34, -44}, {30, -46}, {26, -48}, {22, -50}, {16, -52}, {10, -54}, {4, -54}, {0, -54}, {-8, -54}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 95, 191}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-54, 54}, {56, -54}}, lineColor = {0, 0, 0}, lineThickness = 0.5), Rectangle(extent = {{-4, -52}, {4, -74}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 95, 191}, fillPattern = FillPattern.Solid), Rectangle(extent = {{4, -74}, {88, -68}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 95, 191}, fillPattern = FillPattern.Solid), Line(points = {{-4, -54}, {-4, -74}, {88, -74}}, thickness = 0.5), Line(points = {{4, -54}, {4, -68}, {88, -68}}, thickness = 0.5)}),
      Documentation(info = "This component is used to feed an ASM1 wwtp model with flow data from measurement
  when e.g. concentration is measured after the primary clarifier.
  
  The dimension of InPort is 1.
  
    1 volumeflowrate Q of incoming wastewater [m3/d]"),
      Icon(graphics = {Ellipse(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, extent = {{-54, 54}, {56, -54}}, endAngle = 360), Polygon(lineColor = {0, 0, 255}, fillColor = {0, 95, 191}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-8, -54}, {-14, -52}, {-24, -48}, {-32, -44}, {-36, -40}, {-42, -34}, {-48, -26}, {-50, -20}, {52, -20}, {50, -26}, {46, -32}, {42, -36}, {38, -40}, {34, -44}, {30, -46}, {26, -48}, {22, -50}, {16, -52}, {10, -54}, {4, -54}, {0, -54}, {-8, -54}}), Ellipse(lineThickness = 0.5, extent = {{-54, 54}, {56, -54}}, endAngle = 360), Polygon(origin = {-31, 63}, fillColor = {186, 189, 182}, fillPattern = FillPattern.Solid, points = {{-37, 9}, {37, 9}, {37, -9}, {25, -9}, {25, -3}, {-37, -3}, {-37, 9}})}, coordinateSystem(initialScale = 0.1)));
  end WaterSink;



  annotation(
    uses(Modelica(version = "3.2.2")));
end PumWasteWaterPlant;