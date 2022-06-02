package DieselEngineLibraryTest
  package RadTest
    model HeatExchangerCell1
      import SI = Modelica.SIunits;
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium waterMedium = DieselEngineLibrary.Media.LLC() "Cooling water" annotation(
        choicesAllMatching = true);
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium windMedium = DieselEngineLibrary.Media.Air_75degC() "Cooling wind" annotation(
        choicesAllMatching = true);
      parameter Modelica.SIunits.Length Lrws = 0.525 "wide of sample radiator";
      parameter Modelica.SIunits.Length Lrw = 0.525 "wide of radiator";
      //water
      parameter Modelica.SIunits.VolumeFlowRate V_flow_water = 0 annotation(
        Dialog(group = "initialization"));
      parameter Boolean fixed_V_flow_water = true annotation(
        Dialog(group = "initialization"));
      //  parameter SI.Temperature T_water = 373 annotation(
      //    Dialog(group = "initialization"));
      parameter Modelica.SIunits.Temperature T_water = 373 annotation(
        Dialog(group = "initialization"));
      parameter Boolean fixed_T_water = false annotation(
        Dialog(group = "initialization"));
      parameter Modelica.SIunits.Mass m_water = 0;
      //wind
      parameter Modelica.SIunits.VolumeFlowRate V_flow_wind = 0 annotation(
        Dialog(group = "initialization"));
      parameter Boolean fixed_V_flow_wind = true annotation(
        Dialog(group = "initialization"));
      parameter Modelica.SIunits.Temperature T_wind = 313 annotation(
        Dialog(group = "initialization"));
      parameter Boolean fixed_T_wind = false annotation(
        Dialog(group = "initialization"));
      parameter Modelica.SIunits.Mass m_wind = 0;
      //table data
      parameter String tableName = "Tab1" "Table name on file or in function usertab (see docu)" annotation(
        Dialog(group = "Table data definition"));
      parameter String fileName = "C:\\Work\\2020\\DieselEngineLibrary\\csv\\data_K.txt" "File where matrix is stored" annotation(
        Dialog(group = "Table data definition", loadSelector(filter = "Text files (*.txt);;csv-files (*.csv)", caption = "Open file in which table is present")));
      parameter Modelica.SIunits.Temperature T0_water;
      parameter Modelica.SIunits.Temperature T0_wind;
      //port
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a Ta1(medium = windMedium) annotation(
        Placement(visible = true, transformation(origin = {-142, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b Tw2(medium = waterMedium) annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b Ta2(medium = windMedium) annotation(
        Placement(visible = true, transformation(origin = {140, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //pipe
      //DieselEngineLibrary.Components.Commons.HeatedPipe
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a Tw1(medium = waterMedium) annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Tables.CombiTable2D combiTable2D(fileName = fileName, tableName = tableName, tableOnFile = true) annotation(
        Placement(visible = true, transformation(origin = {-102, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain(k = Lrw / Lrws * 1000) annotation(
        Placement(visible = true, transformation(origin = {-72, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
        Placement(visible = true, transformation(origin = {-36, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Components.HeatedPipe heatedPipe(T0 = T0_wind, T0fixed = true, medium = windMedium) annotation(
        Placement(visible = true, transformation(origin = {-58, 0}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Components.HeatedPipe heatedPipe1(T0 = T0_water, T0fixed = true, medium = waterMedium) annotation(
        Placement(visible = true, transformation(origin = {0, 30}, extent = {{10, 10}, {-10, -10}}, rotation = 90)));
      Modelica.Thermal.FluidHeatFlow.Sensors.VolumeFlowSensor volumeFlowSensor(medium = windMedium) annotation(
        Placement(visible = true, transformation(origin = {-106, 0}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sensors.VolumeFlowSensor volumeFlowSensor1(medium = windMedium) annotation(
        Placement(visible = true, transformation(origin = {0, 72}, extent = {{10, 10}, {-10, -10}}, rotation = 90)));
      //Modelica.SIunits.Temperature TaOut;
      //Modelica.SIunits.Temperature TwOut;
      //algorithm
      //TwOut:=heatedPipe1.T_b;
      //TaOut:=heatedPipe.T_b;
    equation
      connect(combiTable2D.y, gain.u) annotation(
        Line(points = {{-91, 54}, {-84, 54}}, color = {0, 0, 127}));
      connect(gain.y, convection.Gc) annotation(
        Line(points = {{-60, 54}, {-36, 54}, {-36, 40}}, color = {0, 0, 127}));
      connect(heatedPipe.flowPort_b, Ta2) annotation(
        Line(points = {{-48, 0}, {-6, 0}, {-6, 6}, {8, 6}, {8, 0}, {142, 0}, {142, 0}, {140, 0}}, color = {255, 0, 0}));
      connect(heatedPipe.heatPort, convection.solid) annotation(
        Line(points = {{-58, 10}, {-58, 10}, {-58, 30}, {-46, 30}, {-46, 30}}, color = {191, 0, 0}));
      connect(Ta1, volumeFlowSensor.flowPort_a) annotation(
        Line(points = {{-142, 0}, {-118, 0}, {-118, 0}, {-116, 0}}, color = {255, 0, 0}));
      connect(volumeFlowSensor.flowPort_b, heatedPipe.flowPort_a) annotation(
        Line(points = {{-96, 0}, {-68, 0}, {-68, 0}, {-68, 0}}, color = {255, 0, 0}));
      connect(Tw1, volumeFlowSensor1.flowPort_a) annotation(
        Line(points = {{0, 98}, {0, 98}, {0, 82}, {0, 82}}, color = {255, 0, 0}));
      connect(volumeFlowSensor.y, combiTable2D.u1) annotation(
        Line(points = {{-106, 12}, {-128, 12}, {-128, 60}, {-114, 60}, {-114, 60}}, color = {0, 0, 127}));
      connect(combiTable2D.u2, volumeFlowSensor1.y) annotation(
        Line(points = {{-114, 48}, {-124, 48}, {-124, 72}, {-10, 72}, {-10, 72}}, color = {0, 0, 127}));
      connect(volumeFlowSensor1.flowPort_b, heatedPipe1.flowPort_a) annotation(
        Line(points = {{0, 62}, {0, 62}, {0, 40}, {0, 40}}, color = {255, 0, 0}));
      connect(heatedPipe1.heatPort, convection.fluid) annotation(
        Line(points = {{-10, 30}, {-26, 30}, {-26, 30}, {-26, 30}}, color = {191, 0, 0}));
      connect(heatedPipe1.flowPort_b, Tw2) annotation(
        Line(points = {{0, 20}, {0, 20}, {0, -100}, {0, -100}}, color = {255, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-8, -19}, fillColor = {188, 188, 188}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -41}, {60, -61}, {60, 59}, {-60, 79}}), Polygon(origin = {65, -19}, fillColor = {203, 203, 203}, fillPattern = FillPattern.Solid, points = {{-13, 59}, {7, 79}, {7, -41}, {-13, -61}, {-13, -13}, {-13, 59}}), Polygon(origin = {5, 60}, fillColor = {229, 229, 229}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-53, 20}, {67, 0}, {47, -20}, {-73, 0}}), Text(origin = {31, 93}, extent = {{-77, 13}, {11, -27}}, textString = "cooling water"), Text(origin = {-37, -11}, extent = {{-77, 13}, {11, -27}}, textString = "cooling wind")}, coordinateSystem(initialScale = 0.1)),
        Diagram(coordinateSystem(extent = {{-140, -100}, {140, 100}})));
    end HeatExchangerCell1;

    model CellTest2
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell2(T0_water = 373.15, T0_wind = 293.15) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0 = 373.15, constantVolumeFlow = 0.001, m = 0.1, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0 = 293.15, constantVolumeFlow = 0.55, m = 0.1, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(volumeFlow.flowPort_b, heatExchangerCell2.Tw1) annotation(
        Line(points = {{-32, 58}, {-2, 58}, {-2, 10}, {0, 10}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerCell2.Ta1) annotation(
        Line(points = {{-36, 0}, {-12, 0}, {-12, 0}, {-10, 0}}, color = {255, 0, 0}));
      connect(heatExchangerCell2.Ta2, ambient3.flowPort) annotation(
        Line(points = {{10, 0}, {54, 0}, {54, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(ambient2.flowPort, heatExchangerCell2.Tw2) annotation(
        Line(points = {{2, -62}, {0, -62}, {0, -10}, {0, -10}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst -d=bltdump -d=bltdump ",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
    end CellTest2;

    model HeatExchangerBasic1_2
      //parameter
      parameter Integer Nj(min = 1) = 21 "高さ方向の分割数";
      parameter Integer Nk(min = 1) = 18 "奥行方向の分割数";
      parameter Real Lrw(min = 0) = 525 "幅";
      parameter Real Lrh(min = 0) = 525 "高さ";
      parameter Real Lrd(min = 0) = 36 "奥行";
      //port
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a wind_in annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b water_out annotation(
        Placement(visible = true, transformation(origin = {-2, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b wind_out annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //  each V_flow_water = 1e-3/Nk,
      //  each V_flow_wind = 0.55/Nj,
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a flowPort_a annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation

    protected
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-14, -17}, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -57}, {60, -79}, {60, 57}, {-60, 79}}), Polygon(origin = {59, -17}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, points = {{-13, 57}, {13, 79}, {13, -61}, {-13, -79}, {-13, -13}, {-13, 57}}), Polygon(origin = {-1, 62}, fillColor = {218, 218, 218}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-47, 22}, {73, 0}, {47, -22}, {-73, 0}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, 8.275}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.3, -37.2}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.275, -61.5}, points = {{-60, 11}, {60, -11}}), Line(origin = {59, 31}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.025, 8.625}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.85, -14.475}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.875, -36.85}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.075, -61.5}, points = {{-13, -11}, {13, 11}}), Line(origin = {6.825, 69}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-0.85, 63.425}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-8, 57.35}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {54, -21.775}, points = {{0, 69}, {0, -69}}), Line(origin = {60.675, -17.375}, points = {{0, 69}, {0, -69}}), Line(origin = {67.7, -11.75}, points = {{0, 69}, {0, -69}})}, coordinateSystem(initialScale = 0.1)),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump ",
        __OpenModelica_simulationFlags(iim = "symbolic", lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end HeatExchangerBasic1_2;

    model RadTest1_2
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0 = 373.15, constantVolumeFlow = 0.001, m = 0.1, medium = llc, flowPort_a(m_flow(start = 0.001 * llc.rho, fixed = false)), flowPort_b(m_flow(start = -0.001 * llc.rho, fixed = false))) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //  Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0 = 373.15, constantVolumeFlow = 0.001, m = 0.1, medium = llc) annotation(
      //    Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0 = 293.15, constantVolumeFlow = 0.001, m = 0.1, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      HeatExchangerBasic1_2 heatExchangerBasic1_2 annotation(
        Placement(visible = true, transformation(origin = {-2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerBasic1_2.wind_in) annotation(
        Line(points = {{-36, 0}, {-12, 0}, {-12, 0}, {-12, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic1_2.wind_out, ambient3.flowPort) annotation(
        Line(points = {{8, 0}, {56, 0}, {56, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic1_2.water_out, ambient2.flowPort) annotation(
        Line(points = {{-2, -10}, {-6, -10}, {-6, -62}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
      connect(heatExchangerBasic1_2.flowPort_b, ambient2.flowPort) annotation(
        Line(points = {{4, -10}, {2, -10}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
      connect(heatExchangerBasic1_2.flowPort_a, volumeFlow.flowPort_b) annotation(
        Line(points = {{-6, 10}, {-6, 10}, {-6, 58}, {-32, 58}, {-32, 58}}, color = {255, 0, 0}));
      connect(heatExchangerBasic1_2.flowPort_a1, volumeFlow.flowPort_b) annotation(
        Line(points = {{4, 10}, {4, 10}, {4, 58}, {-32, 58}, {-32, 58}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst -d=bltdump -d=bltdump ",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
    end RadTest1_2;

    model HeatExchangerBasic2
      //parameter
      parameter Integer Nj(min = 1) = 21 "高さ方向の分割数";
      parameter Integer Nk(min = 1) = 18 "奥行方向の分割数";
      parameter Real Lrw(min = 0) = 525 "幅";
      parameter Real Lrh(min = 0) = 525 "高さ";
      parameter Real Lrd(min = 0) = 36 "奥行";
      //port
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a wind_in annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a water_in annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b water_out annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b wind_out annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //  each V_flow_water = 1e-3/Nk,
      //  each V_flow_wind = 0.55/Nj,
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell(m_water = 0.1, m_wind = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell1 annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell2 annotation(
        Placement(visible = true, transformation(origin = {60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell3 annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell4(m_water = 0.1, m_wind = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell5 annotation(
        Placement(visible = true, transformation(origin = {60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(wind_in, heatExchangerCell.Tw1) annotation(
        Line(points = {{-100, 0}, {-80, 0}, {-80, 40}, {-70, 40}}, color = {255, 0, 0}));
      connect(wind_in, heatExchangerCell4.Tw1) annotation(
        Line(points = {{-100, 0}, {-80, 0}, {-80, -38}, {-70, -38}, {-70, -40}}, color = {255, 0, 0}));
      connect(heatExchangerCell5.Tw2, wind_out) annotation(
        Line(points = {{70, -40}, {80, -40}, {80, 0}, {100, 0}, {100, 0}}, color = {255, 0, 0}));
      connect(heatExchangerCell2.Tw2, wind_out) annotation(
        Line(points = {{70, 40}, {80, 40}, {80, 0}, {100, 0}, {100, 0}}, color = {255, 0, 0}));
      connect(water_in, heatExchangerCell1.Ta1) annotation(
        Line(points = {{0, 98}, {0, 98}, {0, 50}, {0, 50}}, color = {255, 0, 0}));
      connect(water_in, heatExchangerCell.Ta1) annotation(
        Line(points = {{0, 98}, {0, 98}, {0, 60}, {-60, 60}, {-60, 50}, {-60, 50}}, color = {255, 0, 0}));
      connect(heatExchangerCell.Tw2, heatExchangerCell1.Tw1) annotation(
        Line(points = {{-50, 40}, {-10, 40}, {-10, 40}, {-10, 40}}, color = {255, 0, 0}));
      connect(heatExchangerCell1.Ta2, heatExchangerCell3.Ta1) annotation(
        Line(points = {{0, 30}, {0, 30}, {0, -30}, {0, -30}}, color = {255, 0, 0}));
      connect(heatExchangerCell.Ta2, heatExchangerCell4.Ta1) annotation(
        Line(points = {{-60, 30}, {-60, 30}, {-60, -30}, {-60, -30}}, color = {255, 0, 0}));
      connect(heatExchangerCell4.Tw2, heatExchangerCell3.Tw1) annotation(
        Line(points = {{-50, -40}, {-10, -40}, {-10, -40}, {-10, -40}}, color = {255, 0, 0}));
      connect(heatExchangerCell4.Ta2, water_out) annotation(
        Line(points = {{-60, -50}, {-60, -50}, {-60, -72}, {0, -72}, {0, -100}, {0, -100}}, color = {255, 0, 0}));
      connect(heatExchangerCell3.Ta2, water_out) annotation(
        Line(points = {{0, -50}, {0, -50}, {0, -100}, {0, -100}}, color = {255, 0, 0}));
      connect(heatExchangerCell5.Ta2, water_out) annotation(
        Line(points = {{60, -50}, {60, -50}, {60, -72}, {0, -72}, {0, -100}, {0, -100}}, color = {255, 0, 0}));
      connect(heatExchangerCell3.Tw2, heatExchangerCell5.Tw1) annotation(
        Line(points = {{10, -40}, {50, -40}, {50, -40}, {50, -40}}, color = {255, 0, 0}));
      connect(heatExchangerCell1.Tw2, heatExchangerCell2.Tw1) annotation(
        Line(points = {{10, 40}, {50, 40}, {50, 40}, {50, 40}}, color = {255, 0, 0}));
      connect(heatExchangerCell2.Ta2, heatExchangerCell5.Ta1) annotation(
        Line(points = {{60, 30}, {60, 30}, {60, -30}, {60, -30}}, color = {255, 0, 0}));
      connect(heatExchangerCell2.Ta1, water_in) annotation(
        Line(points = {{60, 50}, {60, 50}, {60, 60}, {0, 60}, {0, 98}, {0, 98}}, color = {255, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-14, -17}, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -57}, {60, -79}, {60, 57}, {-60, 79}}), Polygon(origin = {59, -17}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, points = {{-13, 57}, {13, 79}, {13, -61}, {-13, -79}, {-13, -13}, {-13, 57}}), Polygon(origin = {-1, 62}, fillColor = {218, 218, 218}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-47, 22}, {73, 0}, {47, -22}, {-73, 0}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, 8.275}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.3, -37.2}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.275, -61.5}, points = {{-60, 11}, {60, -11}}), Line(origin = {59, 31}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.025, 8.625}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.85, -14.475}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.875, -36.85}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.075, -61.5}, points = {{-13, -11}, {13, 11}}), Line(origin = {6.825, 69}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-0.85, 63.425}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-8, 57.35}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {54, -21.775}, points = {{0, 69}, {0, -69}}), Line(origin = {60.675, -17.375}, points = {{0, 69}, {0, -69}}), Line(origin = {67.7, -11.75}, points = {{0, 69}, {0, -69}})}, coordinateSystem(initialScale = 0.1)));
    end HeatExchangerBasic2;

    model RadTest2
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0(displayUnit = "K") = 373.15, constantVolumeFlow = 0.001, flowPort_a(m_flow(fixed = false, start = 0.001 * llc.rho)), medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0(displayUnit = "K") = 293.15, constantVolumeFlow = 0.001, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      HeatExchangerBasic2 heatExchangerBasic2 annotation(
        Placement(visible = true, transformation(origin = {2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerBasic2.wind_in) annotation(
        Line(points = {{-36, 0}, {-8, 0}, {-8, -2}, {-8, -2}}, color = {255, 0, 0}));
      connect(volumeFlow.flowPort_b, heatExchangerBasic2.water_in) annotation(
        Line(points = {{-32, 58}, {0, 58}, {0, 8}, {2, 8}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.wind_out, ambient3.flowPort) annotation(
        Line(points = {{12, -2}, {56, -2}, {56, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.water_out, ambient2.flowPort) annotation(
        Line(points = {{2, -12}, {2, -12}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump ",
        __OpenModelica_simulationFlags(iim = "symbolic", lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end RadTest2;

    model RadTest2_2
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0(displayUnit = "K") = 373.15, constantVolumeFlow = 0.001, flowPort_a(m_flow(fixed = false, start = 0.001 * llc.rho)), medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0(displayUnit = "K") = 293.15, constantVolumeFlow = 0.55, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      HeatExchangerBasic2 heatExchangerBasic2 annotation(
        Placement(visible = true, transformation(origin = {2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerBasic2.wind_in) annotation(
        Line(points = {{-36, 0}, {-8, 0}, {-8, -2}, {-8, -2}}, color = {255, 0, 0}));
      connect(volumeFlow.flowPort_b, heatExchangerBasic2.water_in) annotation(
        Line(points = {{-32, 58}, {0, 58}, {0, 8}, {2, 8}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.wind_out, ambient3.flowPort) annotation(
        Line(points = {{12, -2}, {56, -2}, {56, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.water_out, ambient2.flowPort) annotation(
        Line(points = {{2, -12}, {2, -12}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump ",
        __OpenModelica_simulationFlags(iim = "symbolic", lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end RadTest2_2;

    model HeatExchangerBasic3
      //parameter
      parameter Integer Nj(min = 1) = 21 "高さ方向の分割数";
      parameter Integer Nk(min = 1) = 18 "奥行方向の分割数";
      parameter Real Lrw(min = 0) = 525 "幅";
      parameter Real Lrh(min = 0) = 525 "高さ";
      parameter Real Lrd(min = 0) = 36 "奥行";
      //port
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a wind_in annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a water_in annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b water_out annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b wind_out annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //  each V_flow_water = 1e-3/Nk,
      //  each V_flow_wind = 0.55/Nj,
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell(m_water = 0.1, m_wind = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell1 annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell3 annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 heatExchangerCell4(m_water = 0.1, m_wind = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //if Nk == 1
      //  connect(wind_in, cell[1,1].Tw1);
      //  connect(wind_in, cell[2,1].Tw1);
      //  connect(wind_out, cell[1,1].Tw2);
      //  connect(wind_out, cell[2,1].Tw2);
      //  connect(water_in, cell[1,1].Ta1);
      //  connect(water_out, cell[2,1].Ta2);
      //if Nj == 1
      //  connect(wind_in, cell[1,1].Tw1);
      //  connect(wind_out, cell[1,1].Tw2);
      //  connect(water_in, cell[1,1].Ta1);
      //  connect(water_in, cell[1,2].Ta1);
      //  connect(water_out, cell[1,1].Ta2);
      //  connect(water_out, cell[1,2].Ta2);
      //  connect(cell[1, 1].Ta2, cell[2, 1].Ta1);
      //connect(cell[p, q].Ta2, cell[p + 1, q].Ta1);
      //output
    equation
      connect(wind_in, heatExchangerCell.Ta1) annotation(
        Line(points = {{-100, 0}, {-82, 0}, {-82, 40}, {-70, 40}, {-70, 40}}, color = {255, 0, 0}));
      connect(wind_in, heatExchangerCell4.Ta1) annotation(
        Line(points = {{-100, 0}, {-82, 0}, {-82, -40}, {-70, -40}, {-70, -40}}, color = {255, 0, 0}));
      connect(heatExchangerCell4.Ta2, heatExchangerCell3.Ta1) annotation(
        Line(points = {{-50, -40}, {-10, -40}, {-10, -40}, {-10, -40}}, color = {255, 0, 0}));
      connect(heatExchangerCell4.Tw2, water_out) annotation(
        Line(points = {{-60, -50}, {-60, -50}, {-60, -74}, {-2, -74}, {-2, -100}, {0, -100}}, color = {255, 0, 0}));
      connect(heatExchangerCell3.Tw2, water_out) annotation(
        Line(points = {{0, -50}, {0, -50}, {0, -100}, {0, -100}}, color = {255, 0, 0}));
      connect(heatExchangerCell.Ta2, heatExchangerCell1.Ta1) annotation(
        Line(points = {{-50, 40}, {-12, 40}, {-12, 40}, {-10, 40}}, color = {255, 0, 0}));
      connect(heatExchangerCell.Tw2, heatExchangerCell4.Tw1) annotation(
        Line(points = {{-60, 30}, {-60, 30}, {-60, -30}, {-60, -30}}, color = {255, 0, 0}));
      connect(heatExchangerCell.Tw1, water_in) annotation(
        Line(points = {{-60, 50}, {-60, 50}, {-60, 80}, {0, 80}, {0, 98}, {0, 98}}, color = {255, 0, 0}));
      connect(water_in, heatExchangerCell1.Tw1) annotation(
        Line(points = {{0, 98}, {0, 98}, {0, 50}, {0, 50}}, color = {255, 0, 0}));
      connect(heatExchangerCell1.Tw2, heatExchangerCell3.Tw1) annotation(
        Line(points = {{0, 30}, {0, 30}, {0, -30}, {0, -30}}, color = {255, 0, 0}));
      connect(heatExchangerCell1.Ta2, wind_out) annotation(
        Line(points = {{10, 40}, {40, 40}, {40, 0}, {100, 0}, {100, 0}}, color = {255, 0, 0}));
      connect(heatExchangerCell3.Ta2, wind_out) annotation(
        Line(points = {{10, -40}, {40, -40}, {40, 0}, {100, 0}, {100, 0}}, color = {255, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-14, -17}, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -57}, {60, -79}, {60, 57}, {-60, 79}}), Polygon(origin = {59, -17}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, points = {{-13, 57}, {13, 79}, {13, -61}, {-13, -79}, {-13, -13}, {-13, 57}}), Polygon(origin = {-1, 62}, fillColor = {218, 218, 218}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-47, 22}, {73, 0}, {47, -22}, {-73, 0}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, 8.275}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.3, -37.2}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.275, -61.5}, points = {{-60, 11}, {60, -11}}), Line(origin = {59, 31}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.025, 8.625}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.85, -14.475}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.875, -36.85}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.075, -61.5}, points = {{-13, -11}, {13, 11}}), Line(origin = {6.825, 69}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-0.85, 63.425}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-8, 57.35}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {54, -21.775}, points = {{0, 69}, {0, -69}}), Line(origin = {60.675, -17.375}, points = {{0, 69}, {0, -69}}), Line(origin = {67.7, -11.75}, points = {{0, 69}, {0, -69}})}, coordinateSystem(initialScale = 0.1)));
    end HeatExchangerBasic3;

    model RadTest3
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0 = 373.15, constantVolumeFlow = 0.001, m = 0.1, medium = llc, flowPort_a(m_flow(start = 0.001 * llc.rho, fixed = true))) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0 = 293.15, constantVolumeFlow = 0.001, m = 0.1, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      HeatExchangerBasic3 heatExchangerBasic2 annotation(
        Placement(visible = true, transformation(origin = {2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerBasic2.wind_in) annotation(
        Line(points = {{-36, 0}, {-8, 0}, {-8, -2}, {-8, -2}}, color = {255, 0, 0}));
      connect(volumeFlow.flowPort_b, heatExchangerBasic2.water_in) annotation(
        Line(points = {{-32, 58}, {0, 58}, {0, 8}, {2, 8}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.wind_out, ambient3.flowPort) annotation(
        Line(points = {{12, -2}, {56, -2}, {56, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.water_out, ambient2.flowPort) annotation(
        Line(points = {{2, -12}, {2, -12}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
    end RadTest3;

    model HeatExchangerBasic4
      //parameter
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium waterMedium = DieselEngineLibrary.Media.LLC() "Cooling water" annotation(
        choicesAllMatching = true);
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium windMedium = DieselEngineLibrary.Media.Air_75degC() "Cooling wind" annotation(
        choicesAllMatching = true);
      parameter Integer Nj(min = 1) = 2 "高さ方向の分割数";
      parameter Integer Nk(min = 1) = 1 "奥行方向の分割数";
      parameter Real Lrw(min = 0) = 525 "幅";
      parameter Real Lrh(min = 0) = 525 "高さ";
      parameter Real Lrd(min = 0) = 36 "奥行";
      parameter Modelica.SIunits.Temperature T0_water;
      parameter Modelica.SIunits.Temperature T0_wind;
      //variable
      Modelica.SIunits.Temperature Tin_water;
      Modelica.SIunits.Temperature Tout_water;
      Modelica.SIunits.Temperature Tin_wind;
      Modelica.SIunits.Temperature Tout_wind;
      //port
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a wind_in annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a water_in annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b water_out annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b wind_out annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //  DieselEngineLibrary.Components.Radiator.HeatExchangerCell1 cell[2]  annotation(
      //    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerCell1 cell[Nj, Nk] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    protected
      parameter Integer Mj = max(1, Nj - 1);
      parameter Integer Mk = max(1, Nk - 1);
    equation
      for n in 1:Nj loop
        connect(wind_in, cell[n, 1].Tw1);
        connect(wind_out, cell[n, Nk].Tw2);
      end for;
      for m in 1:Nk loop
        connect(water_in, cell[1, m].Ta1);
        connect(water_out, cell[Nj, m].Ta2);
      end for;
      if Nj > 1 and Nk == 1 then
        for p in 1:Mj loop
          for q in 1:Mk loop
            connect(cell[p, q].Ta2, cell[p + 1, q].Ta1);
          end for;
        end for;
      elseif Nj == 1 and Nk > 1 then
        for p in 1:Mj loop
          for q in 1:Mk loop
            connect(cell[p, q].Tw2, cell[p, q + 1].Tw1);
          end for;
        end for;
      else
        for p in 1:Mj loop
          for q in 1:Mk loop
            connect(cell[p, q].Tw2, cell[p, q + 1].Tw1);
            connect(cell[p, q].Ta2, cell[p + 1, q].Ta1);
          end for;
        end for;
      end if;
      Tin_water = water_in.h / waterMedium.cp;
      Tout_water = water_out.h / waterMedium.cp;
      Tin_wind = wind_in.h / windMedium.cp;
      Tout_wind = wind_out.h / windMedium.cp;
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-14, -17}, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -57}, {60, -79}, {60, 57}, {-60, 79}}), Polygon(origin = {59, -17}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, points = {{-13, 57}, {13, 79}, {13, -61}, {-13, -79}, {-13, -13}, {-13, 57}}), Polygon(origin = {-1, 62}, fillColor = {218, 218, 218}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-47, 22}, {73, 0}, {47, -22}, {-73, 0}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, 8.275}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.3, -37.2}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.275, -61.5}, points = {{-60, 11}, {60, -11}}), Line(origin = {59, 31}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.025, 8.625}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.85, -14.475}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.875, -36.85}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.075, -61.5}, points = {{-13, -11}, {13, 11}}), Line(origin = {6.825, 69}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-0.85, 63.425}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-8, 57.35}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {54, -21.775}, points = {{0, 69}, {0, -69}}), Line(origin = {60.675, -17.375}, points = {{0, 69}, {0, -69}}), Line(origin = {67.7, -11.75}, points = {{0, 69}, {0, -69}})}));
    end HeatExchangerBasic4;

    model RadTest4
      parameter DieselEngineLibrary.Media.LLC llc(cp = 3320) annotation(
        Placement(visible = true, transformation(origin = {-46, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter DieselEngineLibrary.Media.Air_75degC air annotation(
        Placement(visible = true, transformation(origin = {-72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(T0 = 373.15, constantVolumeFlow = 0.001, m = 0.1, medium = llc, flowPort_a(m_flow(start = 0.001 * llc.rho, fixed = true))) annotation(
        Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-78, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient1(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-82, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow1(T0 = 293.15, constantVolumeFlow = 0.001, m = 0.1, medium = air) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient3(constantAmbientPressure = 100000, constantAmbientTemperature = 293.15, medium = air) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Sources.Ambient ambient2(constantAmbientPressure = 100000, constantAmbientTemperature = 373.15, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibraryTest.HeatExchangerBasic4 heatExchangerBasic2(Nk = 2) annotation(
        Placement(visible = true, transformation(origin = {2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ambient.flowPort, volumeFlow.flowPort_a) annotation(
        Line(points = {{-68, 58}, {-52, 58}, {-52, 58}, {-52, 58}}, color = {255, 0, 0}));
      connect(ambient1.flowPort, volumeFlow1.flowPort_a) annotation(
        Line(points = {{-72, 0}, {-56, 0}, {-56, 0}, {-56, 0}}, color = {255, 0, 0}));
      connect(volumeFlow1.flowPort_b, heatExchangerBasic2.wind_in) annotation(
        Line(points = {{-36, 0}, {-8, 0}, {-8, -2}, {-8, -2}}, color = {255, 0, 0}));
      connect(volumeFlow.flowPort_b, heatExchangerBasic2.water_in) annotation(
        Line(points = {{-32, 58}, {0, 58}, {0, 8}, {2, 8}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.wind_out, ambient3.flowPort) annotation(
        Line(points = {{12, -2}, {56, -2}, {56, 0}, {56, 0}}, color = {255, 0, 0}));
      connect(heatExchangerBasic2.water_out, ambient2.flowPort) annotation(
        Line(points = {{2, -12}, {2, -12}, {2, -62}, {2, -62}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
    end RadTest4;
    annotation(
      uses(Modelica(version = "3.2.3")));
  end RadTest;

  package ConfluenceTest
    model PipeConfluence2_4
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {26, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow1(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {22, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {8, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.FixedVolumeFlow1 fixedVolumeFlow1(T_a(displayUnit = "K") = 293.15, V_flow = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-98, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.Ambient2 ambient21 annotation(
        Placement(visible = true, transformation(origin = {100, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow1(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {-4, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ConfluenceOneToMulti confluenceOneToMulti(m_flowIni = -0.05) annotation(
        Placement(visible = true, transformation(origin = {-58, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      ConfluenceMultiToOne confluenceMultiToOne annotation(
        Placement(visible = true, transformation(origin = {64, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      //,a2.m_flow(start=-0.05,fixed=true),a1.m_flow(start=-0.05,fixed=true)
    equation
      connect(fixedHeatFlow.port, heatedPipeFixedVolumeFlow.port_a) annotation(
        Line(points = {{18, 18}, {18, 10}, {26, 10}, {26, 2}}, color = {191, 0, 0}));
      connect(fixedHeatFlow1.port, heatedPipeFixedVolumeFlow1.port_a) annotation(
        Line(points = {{6, -38}, {22, -38}, {22, -62}}, color = {191, 0, 0}));
      connect(fixedVolumeFlow1.flowPort_a, confluenceOneToMulti.b) annotation(
        Line(points = {{-88, -42}, {-68, -42}, {-68, -40}}, color = {255, 0, 0}));
      connect(confluenceOneToMulti.a2, heatedPipeFixedVolumeFlow.flowPort_a) annotation(
        Line(points = {{-48, -34}, {-26, -34}, {-26, -8}, {16, -8}}, color = {255, 0, 0}));
      connect(confluenceOneToMulti.a1, heatedPipeFixedVolumeFlow1.flowPort_a) annotation(
        Line(points = {{-48, -46}, {-28, -46}, {-28, -72}, {12, -72}}, color = {255, 0, 0}));
      connect(heatedPipeFixedVolumeFlow.flowPort_b, confluenceMultiToOne.a2) annotation(
        Line(points = {{36, -8}, {42, -8}, {42, -38}, {54, -38}, {54, -38}}, color = {255, 0, 0}));
      connect(heatedPipeFixedVolumeFlow1.flowPort_b, confluenceMultiToOne.a1) annotation(
        Line(points = {{32, -72}, {44, -72}, {44, -50}, {54, -50}, {54, -50}}, color = {255, 0, 0}));
      connect(confluenceMultiToOne.b, ambient21.flowPort) annotation(
        Line(points = {{74, -44}, {88, -44}, {88, -44}, {90, -44}}, color = {255, 0, 0}));
    protected
      annotation(
        uses(Modelica(version = "3.2.3")));
    end PipeConfluence2_4;

    model PipeConfluence2_3
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {26, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow1(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {22, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {8, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.FixedVolumeFlow1 fixedVolumeFlow1(T_a(displayUnit = "K") = 293.15, V_flow = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-98, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.Ambient2 ambient21 annotation(
        Placement(visible = true, transformation(origin = {54, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.Ambient2 ambient2 annotation(
        Placement(visible = true, transformation(origin = {56, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow1(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {-4, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ConfluenceOneToMulti confluenceOneToMulti(m_flowIni = -0.05) annotation(
        Placement(visible = true, transformation(origin = {-58, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      //,a2.m_flow(start=-0.05,fixed=true),a1.m_flow(start=-0.05,fixed=true)
    equation
      connect(fixedHeatFlow.port, heatedPipeFixedVolumeFlow.port_a) annotation(
        Line(points = {{18, 18}, {18, 10}, {26, 10}, {26, 2}}, color = {191, 0, 0}));
      connect(heatedPipeFixedVolumeFlow1.flowPort_b, ambient21.flowPort) annotation(
        Line(points = {{32, -72}, {44, -72}}, color = {255, 0, 0}));
      connect(heatedPipeFixedVolumeFlow.flowPort_b, ambient2.flowPort) annotation(
        Line(points = {{36, -8}, {46, -8}}, color = {255, 0, 0}));
      connect(fixedHeatFlow1.port, heatedPipeFixedVolumeFlow1.port_a) annotation(
        Line(points = {{6, -38}, {22, -38}, {22, -62}}, color = {191, 0, 0}));
      connect(fixedVolumeFlow1.flowPort_a, confluenceOneToMulti.b) annotation(
        Line(points = {{-88, -42}, {-68, -42}, {-68, -40}}, color = {255, 0, 0}));
      connect(confluenceOneToMulti.a2, heatedPipeFixedVolumeFlow.flowPort_a) annotation(
        Line(points = {{-48, -34}, {-26, -34}, {-26, -8}, {16, -8}}, color = {255, 0, 0}));
      connect(confluenceOneToMulti.a1, heatedPipeFixedVolumeFlow1.flowPort_a) annotation(
        Line(points = {{-48, -46}, {-28, -46}, {-28, -72}, {12, -72}}, color = {255, 0, 0}));
    protected
      annotation(
        uses(Modelica(version = "3.2.3")));
    end PipeConfluence2_3;

    model ConfluenceOneToMulti
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a a1(m_flow(start = m_flowIni, fixed = true)) annotation(
        Placement(visible = true, transformation(origin = {-60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a a2(m_flow(start = m_flowIni, fixed = true)) annotation(
        Placement(visible = true, transformation(origin = {60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b b annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter Integer n = 2;
      parameter Real m_flowIni;
      //initial algorithm
      //  a2.m_flow:=a2m_flowIni;
      //  a1.m_flow:=a2m_flowIni;
    algorithm
      a1.H_flow := -b.H_flow / n;
      a2.H_flow := -b.H_flow / n;
      a1.m_flow := -(b.m_flow + 1e-10) / n;
      a2.m_flow := -(b.m_flow + 1e-10) / n;
      b.h := b.H_flow / b.m_flow;
      b.p := a1.p;
      b.p := a2.p;
    equation
//  a1.H_flow=-b.H_flow/n;
//  a2.H_flow=-b.H_flow/n;
//  a1.m_flow=-b.m_flow/n;
//  a2.m_flow=-b.m_flow/n;
//  b.h=b.H_flow/b.m_flow;
//  a2.p = b.p;
//  a1.p = b.p;
      annotation(
        Diagram,
        uses(Modelica(version = "3.2.3")),
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Line(origin = {-30.8867, 10.2013}, points = {{30.8867, 69.7987}, {-29.1133, -70.2013}, {-7.11332, -60.2013}, {-29.1133, -70.2013}, {-31.1133, -48.2013}}), Line(origin = {31.8367, 10.2013}, points = {{-31.8367, 69.7987}, {28.1633, -70.2013}, {32.1633, -46.2013}, {28.1633, -70.2013}, {6.16334, -60.2013}})}));
    end ConfluenceOneToMulti;

    model ConfluenceMultiToOne
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a a1 annotation(
        Placement(visible = true, transformation(origin = {-60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_a a2 annotation(
        Placement(visible = true, transformation(origin = {60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {58, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.FluidHeatFlow.Interfaces.FlowPort_b b annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      a1.H_flow + a2.H_flow + b.H_flow = 0.0;
//  a1.h=a1.H_flow/a1.m_flow;
//  a2.h=a2.H_flow/a2.m_flow;
      b.h = b.H_flow / b.m_flow;
      a1.m_flow + a2.m_flow + b.m_flow = 0.0;
      a2.p = b.p;
      a1.p = b.p;
      annotation(
        Diagram,
        uses(Modelica(version = "3.2.3")),
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Line(origin = {-39.4514, 9.80155}, points = {{-31.8367, 69.7987}, {28.1633, -70.2013}, {32.1633, -46.2013}, {28.1633, -70.2013}, {6.16334, -60.2013}}), Line(origin = {40.6679, 9.80155}, points = {{30.8867, 69.7987}, {-29.1133, -70.2013}, {-7.11332, -60.2013}, {-29.1133, -70.2013}, {-31.1133, -48.2013}})}));
    end ConfluenceMultiToOne;

    model PipeConfluence2_2
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {-34, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.FixedVolumeFlow1 fixedVolumeFlow11(T_a = 293.15, V_flow = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-88, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.HeatedPipeNoReverse heatedPipeFixedVolumeFlow1(T0 = 293.15) annotation(
        Placement(visible = true, transformation(origin = {-42, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {-90, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.FixedVolumeFlow1 fixedVolumeFlow1(T_a = 293.15, V_flow = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-88, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow1(Q_flow = 10) annotation(
        Placement(visible = true, transformation(origin = {-86, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Components.Radiator.Ambient2 ambient21 annotation(
        Placement(visible = true, transformation(origin = {46, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ConfluenceMultiToOne confluenceMultiToOne annotation(
        Placement(visible = true, transformation(origin = {6, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(fixedVolumeFlow11.flowPort_a, heatedPipeFixedVolumeFlow1.flowPort_a) annotation(
        Line(points = {{-78, -72}, {-52, -72}}, color = {255, 0, 0}));
      connect(fixedVolumeFlow1.flowPort_a, heatedPipeFixedVolumeFlow.flowPort_a) annotation(
        Line(points = {{-78, 0}, {-64, 0}, {-64, -8}, {-44, -8}}, color = {255, 0, 0}));
      connect(fixedHeatFlow.port, heatedPipeFixedVolumeFlow.port_a) annotation(
        Line(points = {{-80, 32}, {-44, 32}, {-44, 2}, {-34, 2}}, color = {191, 0, 0}));
      connect(fixedHeatFlow1.port, heatedPipeFixedVolumeFlow1.port_a) annotation(
        Line(points = {{-76, -40}, {-42, -40}, {-42, -62}, {-42, -62}}, color = {191, 0, 0}));
      connect(confluenceMultiToOne.b, ambient21.flowPort) annotation(
        Line(points = {{16, -38}, {24, -38}, {24, -46}, {36, -46}}, color = {255, 0, 0}));
      connect(heatedPipeFixedVolumeFlow.flowPort_b, confluenceMultiToOne.a2) annotation(
        Line(points = {{-24, -8}, {-16, -8}, {-16, -32}, {-4, -32}}, color = {255, 0, 0}));
      connect(heatedPipeFixedVolumeFlow1.flowPort_b, confluenceMultiToOne.a1) annotation(
        Line(points = {{-32, -72}, {-16, -72}, {-16, -44}, {-4, -44}}, color = {255, 0, 0}));
    protected
      annotation(
        uses(Modelica(version = "3.2.3")));
    end PipeConfluence2_2;
  end ConfluenceTest;

  package RadTest2_CustomPipe
  end RadTest2_CustomPipe;
end DieselEngineLibraryTest;
