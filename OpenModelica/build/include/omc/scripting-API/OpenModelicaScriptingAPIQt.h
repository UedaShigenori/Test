#include <QtCore>
#include "OpenModelicaScriptingAPI.h"

class OMCInterface : public QObject
{
  Q_OBJECT
public:
  threadData_t *threadData;
  OMCInterface(threadData_t *td);
  QString oms_getVersion();
  modelica_integer oms_terminate(QString cref);
  modelica_integer oms_stepUntil(QString cref, modelica_real stopTime);
  modelica_integer oms_simulate(QString cref);
  modelica_integer oms_setWorkingDirectory(QString newWorkingDir);
  modelica_integer oms_setVariableStepSize(QString cref, modelica_real initialStepSize, modelica_real minimumStepSize, modelica_real maximumStepSize);
  modelica_integer oms_setTolerance(QString cref, modelica_real absoluteTolerance, modelica_real relativeTolerance);
  modelica_integer oms_setTLMSocketData(QString cref, QString address, modelica_integer managerPort, modelica_integer monitorPort);
  modelica_integer oms_setTLMPositionAndOrientation(QString cref, modelica_real x1, modelica_real x2, modelica_real x3, modelica_real A11, modelica_real A12, modelica_real A13, modelica_real A21, modelica_real A22, modelica_real A23, modelica_real A31, modelica_real A32, modelica_real A33);
  modelica_integer oms_setTempDirectory(QString newTempDir);
  modelica_integer oms_setStopTime(QString cref, modelica_real stopTime);
  modelica_integer oms_setStartTime(QString cref, modelica_real startTime);
  modelica_integer oms_setSignalFilter(QString cref, QString regex);
  modelica_integer oms_setResultFile(QString cref, QString filename, modelica_integer bufferSize);
  modelica_integer oms_setRealInputDerivative(QString cref, modelica_real value);
  modelica_integer oms_setReal(QString cref, modelica_real value);
  modelica_integer oms_setLoggingLevel(modelica_integer logLevel);
  modelica_integer oms_setLoggingInterval(QString cref, modelica_real loggingInterval);
  modelica_integer oms_setLogFile(QString filename);
  modelica_integer oms_setInteger(QString cref, modelica_integer value);
  modelica_integer oms_setFixedStepSize(QString cref, modelica_real stepSize);
  modelica_integer oms_setCommandLineOption(QString cmd);
  modelica_integer oms_setBoolean(QString cref, modelica_boolean value);
  modelica_integer oms_RunFile(QString filename);
  modelica_integer oms_reset(QString cref);
  modelica_integer oms_rename(QString cref, QString newCref);
  modelica_integer oms_removeSignalsFromResults(QString cref, QString regex);
  typedef struct {
    modelica_integer status;
    QString cref;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append("\"" + cref + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_parseModelName_res;
  oms_parseModelName_res oms_parseModelName(QString contents);
  modelica_integer oms_newModel(QString cref);
  modelica_integer oms_loadSnapshot(QString cref, QString snapshot);
  typedef struct {
    modelica_integer status;
    QString contents;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append("\"" + contents + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_listUnconnectedConnectors_res;
  oms_listUnconnectedConnectors_res oms_listUnconnectedConnectors(QString cref);
  typedef struct {
    modelica_integer status;
    QString contents;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append("\"" + contents + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_list_res;
  oms_list_res oms_list(QString cref);
  modelica_integer oms_instantiate(QString cref);
  modelica_integer oms_initialize(QString cref);
  typedef struct {
    modelica_integer status;
    QString cref;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append("\"" + cref + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_importFile_res;
  oms_importFile_res oms_importFile(QString filename);
  typedef struct {
    modelica_integer status;
    modelica_real initialStepSize;
    modelica_real minimumStepSize;
    modelica_real maximumStepSize;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(initialStepSize));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(minimumStepSize));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(maximumStepSize));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getVariableStepSize_res;
  oms_getVariableStepSize_res oms_getVariableStepSize(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_real absoluteTolerance;
    modelica_real relativeTolerance;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(absoluteTolerance));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(relativeTolerance));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getTolerance_res;
  oms_getTolerance_res oms_getTolerance(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_integer type_;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(type_));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getSystemType_res;
  oms_getSystemType_res oms_getSystemType(QString cref);
  typedef struct {
    modelica_integer status;
    QString path;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append("\"" + path + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getSubModelPath_res;
  oms_getSubModelPath_res oms_getSubModelPath(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_real stopTime;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(stopTime));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getStopTime_res;
  oms_getStopTime_res oms_getStopTime(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_real startTime;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(startTime));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getStartTime_res;
  oms_getStartTime_res oms_getStartTime(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_integer solver;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(solver));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getSolver_res;
  oms_getSolver_res oms_getSolver(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_real value;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(value));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getReal_res;
  oms_getReal_res oms_getReal(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_integer modelState;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(modelState));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getModelState_res;
  oms_getModelState_res oms_getModelState(QString cref);
  modelica_integer oms_getInteger(QString cref, modelica_integer value);
  typedef struct {
    modelica_integer status;
    modelica_real stepSize;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(stepSize));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getFixedStepSize_res;
  oms_getFixedStepSize_res oms_getFixedStepSize(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_boolean value;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(value ? "true" : "false");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_getBoolean_res;
  oms_getBoolean_res oms_getBoolean(QString cref);
  typedef struct {
    modelica_integer status;
    modelica_integer kind;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(status));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(kind));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } oms_extractFMIKind_res;
  oms_extractFMIKind_res oms_extractFMIKind(QString filename);
  modelica_integer oms_exportDependencyGraphs(QString cref, QString initialization, QString simulation);
  modelica_integer oms_export(QString cref, QString filename);
  modelica_integer oms_deleteConnectorFromTLMBus(QString busCref, QString connectorCref);
  modelica_integer oms_deleteConnectorFromBus(QString busCref, QString connectorCref);
  modelica_integer oms_deleteConnection(QString crefA, QString crefB);
  modelica_integer oms_delete(QString cref);
  modelica_integer oms_copySystem(QString source, QString target);
  modelica_integer oms_compareSimulationResults(QString filenameA, QString filenameB, QString var, modelica_real relTol, modelica_real absTol);
  modelica_integer oms_cancelSimulation_asynchronous(QString cref);
  modelica_integer oms_addTLMConnection(QString crefA, QString crefB, modelica_real delay, modelica_real alpha, modelica_real linearimpedance, modelica_real angularimpedance);
  modelica_integer oms_addTimeIndicator(QString signal);
  modelica_integer oms_addSubModel(QString cref, QString fmuPath);
  modelica_integer oms_addStaticValueIndicator(QString signal, modelica_real lower, modelica_real upper, modelica_real stepSize);
  modelica_integer oms_addSignalsToResults(QString cref, QString regex);
  modelica_integer oms_addExternalModel(QString cref, QString path, QString startscript);
  modelica_integer oms_addEventIndicator(QString signal);
  modelica_integer oms_addDynamicValueIndicator(QString signal, QString lower, QString upper, modelica_real stepSize);
  modelica_integer oms_addConnectorToTLMBus(QString busCref, QString connectorCref, QString type_);
  modelica_integer oms_addConnectorToBus(QString busCref, QString connectorCref);
  modelica_integer oms_addConnection(QString crefA, QString crefB);
  modelica_integer oms_addBus(QString cref);
  modelica_integer unloadOMSimulator();
  modelica_integer loadOMSimulator();
  typedef struct {
    modelica_boolean success;
    QString moFile;
    QString qtFile;
    QString qtHeader;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(success ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("\"" + moFile + "\"");
      resultBuffer.append(",");
      resultBuffer.append("\"" + qtFile + "\"");
      resultBuffer.append(",");
      resultBuffer.append("\"" + qtHeader + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } generateScriptingAPI_res;
  generateScriptingAPI_res generateScriptingAPI(QString cl, QString name);
  modelica_boolean deleteInitialState(QString cl, QString state);
  QList<QList<QString > > getInitialStates(QString cl);
  modelica_boolean deleteTransition(QString cl, QString from, QString to, QString condition, modelica_boolean immediate, modelica_boolean reset, modelica_boolean synchronize, modelica_integer priority);
  QList<QList<QString > > getTransitions(QString cl);
  typedef struct {
    QString restriction;
    QString comment;
    modelica_boolean partialPrefix;
    modelica_boolean finalPrefix;
    modelica_boolean encapsulatedPrefix;
    QString fileName;
    modelica_boolean fileReadOnly;
    modelica_integer lineNumberStart;
    modelica_integer columnNumberStart;
    modelica_integer lineNumberEnd;
    modelica_integer columnNumberEnd;
    QList<QString > dimensions;
    modelica_boolean isProtectedClass;
    modelica_boolean isDocumentationClass;
    QString version;
    QString preferredView;
    modelica_boolean state;
    QString access;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append("\"" + restriction + "\"");
      resultBuffer.append(",");
      resultBuffer.append("\"" + comment + "\"");
      resultBuffer.append(",");
      resultBuffer.append(partialPrefix ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(finalPrefix ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(encapsulatedPrefix ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("\"" + fileName + "\"");
      resultBuffer.append(",");
      resultBuffer.append(fileReadOnly ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(QString::number(lineNumberStart));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(columnNumberStart));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(lineNumberEnd));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(columnNumberEnd));
      resultBuffer.append(",");
      resultBuffer.append("{");
      int dimensions_i = 0;
      foreach(QString dimensions_elt, dimensions) {
        if (dimensions_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + dimensions_elt + "\"");
        dimensions_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(",");
      resultBuffer.append(isProtectedClass ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(isDocumentationClass ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("\"" + version + "\"");
      resultBuffer.append(",");
      resultBuffer.append("\"" + preferredView + "\"");
      resultBuffer.append(",");
      resultBuffer.append(state ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("\"" + access + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getClassInformation_res;
  getClassInformation_res getClassInformation(QString cl);
  QList<QString > sortStrings(QList<QString > arr);
  modelica_boolean checkInterfaceOfPackages(QString cl, QList<QList<QString > > dependencyMatrix);
  modelica_boolean GC_set_max_heap_size(modelica_integer size);
  modelica_boolean GC_expand_hp(modelica_integer size);
  void GC_gcollect_and_unmap();
  modelica_real getMemorySize();
  void threadWorkFailed();
  void exit(modelica_integer status);
  QList<modelica_boolean > runScriptParallel(QList<QString > scripts, modelica_integer numThreads, modelica_boolean useThreads);
  modelica_integer numProcessors();
  void generateEntryPoint(QString fileName, QString entryPoint, QString url);
  QString getDerivedClassModifierValue(QString className, QString modifierName);
  QList<QString > getDerivedClassModifierNames(QString className);
  QList<QList<QString > > getUses(QString pack);
  QList<QString > getAvailableLibraries();
  QList<QString > searchClassNames(QString searchText, modelica_boolean findInText);
  modelica_boolean extendsFrom(QString className, QString baseClassName);
  modelica_boolean getBooleanClassAnnotation(QString className, QString annotationName);
  modelica_boolean classAnnotationExists(QString className, QString annotationName);
  QString getAnnotationModifierValue(QString name, QString vendorannotation, QString modifiername);
  QList<QString > getAnnotationNamedModifiers(QString name, QString vendorannotation);
  typedef struct {
    modelica_real startTime;
    modelica_real stopTime;
    modelica_real tolerance;
    modelica_integer numberOfIntervals;
    modelica_real interval;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(startTime));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(stopTime));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(tolerance));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(numberOfIntervals));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(interval));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getSimulationOptions_res;
  getSimulationOptions_res getSimulationOptions(QString name, modelica_real defaultStartTime, modelica_real defaultStopTime, modelica_real defaultTolerance, modelica_integer defaultNumberOfIntervals, modelica_real defaultInterval);
  modelica_boolean isExperiment(QString name);
  QList<QString > getInheritedClasses(QString name);
  QString getBuiltinType(QString cl);
  modelica_boolean isProtectedClass(QString cl, QString c2);
  modelica_boolean isOperatorFunction(QString cl);
  modelica_boolean isOperatorRecord(QString cl);
  modelica_boolean isOperator(QString cl);
  modelica_boolean isEnumeration(QString cl);
  modelica_boolean isOptimization(QString cl);
  modelica_boolean isConnector(QString cl);
  modelica_boolean isModel(QString cl);
  modelica_boolean isPartial(QString cl);
  modelica_boolean isFunction(QString cl);
  modelica_boolean isBlock(QString cl);
  modelica_boolean isRecord(QString cl);
  modelica_boolean isClass(QString cl);
  modelica_boolean isPackage(QString cl);
  modelica_boolean isType(QString cl);
  QString getClassRestriction(QString cl);
  QString basename(QString path);
  QString dirname(QString path);
  QString getClassComment(QString cl);
  QList<QString > typeNameStrings(QString cl);
  QString typeNameString(QString cl);
  QString stringTypeName(QString str);
  typedef struct {
    modelica_real timeStamp;
    QString timeStampAsString;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(timeStamp));
      resultBuffer.append(",");
      resultBuffer.append("\"" + timeStampAsString + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getTimeStamp_res;
  getTimeStamp_res getTimeStamp(QString cl);
  modelica_boolean setDocumentationAnnotation(QString class_, QString info, QString revisions);
  QList<QString > getDocumentationAnnotation(QString cl);
  QString iconv(QString string, QString from, QString to);
  QList<QString > getNthImport(QString class_, modelica_integer index);
  modelica_integer getImportCount(QString class_);
  QString getNthAnnotationString(QString class_, modelica_integer index);
  modelica_integer getAnnotationCount(QString class_);
  QString getNthInitialEquationItem(QString class_, modelica_integer index);
  modelica_integer getInitialEquationItemsCount(QString class_);
  QString getNthEquationItem(QString class_, modelica_integer index);
  modelica_integer getEquationItemsCount(QString class_);
  QString getNthInitialEquation(QString class_, modelica_integer index);
  modelica_integer getInitialEquationCount(QString class_);
  QString getNthEquation(QString class_, modelica_integer index);
  modelica_integer getEquationCount(QString class_);
  QString getNthInitialAlgorithmItem(QString class_, modelica_integer index);
  modelica_integer getInitialAlgorithmItemsCount(QString class_);
  QString getNthAlgorithmItem(QString class_, modelica_integer index);
  modelica_integer getAlgorithmItemsCount(QString class_);
  QString getNthInitialAlgorithm(QString class_, modelica_integer index);
  modelica_integer getInitialAlgorithmCount(QString class_);
  QString getNthAlgorithm(QString class_, modelica_integer index);
  modelica_integer getAlgorithmCount(QString class_);
  QList<QString > getNthConnection(QString className, modelica_integer index);
  modelica_integer getConnectionCount(QString className);
  modelica_boolean updateConnectionNames(QString className, QString from, QString to, QString fromNew, QString toNew);
  modelica_boolean removeExtendsModifiers(QString className, QString baseClassName, modelica_boolean keepRedeclares);
  modelica_boolean removeComponentModifiers(QString class_, QString componentName, modelica_boolean keepRedeclares);
  QList<QString > getInstantiatedParametersAndValues(QString cls);
  QString getComponentModifierValues(QString class_, QString modifier);
  QString getComponentModifierValue(QString class_, QString modifier);
  QList<QString > getComponentModifierNames(QString class_, QString componentName);
  QString getParameterValue(QString class_, QString parameterName);
  QList<QString > getParameterNames(QString class_);
  modelica_boolean closeSimulationResultFile();
  QList<QString > checkCodeGraph(QString graphfile, QString codefile);
  QList<QString > checkTaskGraph(QString filename, QString reffilename);
  QString diffSimulationResultsHtml(QString var, QString actualFile, QString expectedFile, modelica_real relTol, modelica_real relTolDiffMinMax, modelica_real rangeDelta);
  typedef struct {
    modelica_boolean success;
    QList<QString > failVars;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(success ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("{");
      int failVars_i = 0;
      foreach(QString failVars_elt, failVars) {
        if (failVars_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + failVars_elt + "\"");
        failVars_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } diffSimulationResults_res;
  diffSimulationResults_res diffSimulationResults(QString actualFile, QString expectedFile, QString diffPrefix, modelica_real relTol, modelica_real relTolDiffMinMax, modelica_real rangeDelta, QList<QString > vars, modelica_boolean keepEqualResults);
  modelica_real deltaSimulationResults(QString filename, QString reffilename, QString method, QList<QString > vars);
  QList<QString > compareSimulationResults(QString filename, QString reffilename, QString logfilename, modelica_real relTol, modelica_real absTol, QList<QString > vars);
  modelica_boolean filterSimulationResults(QString inFile, QString outFile, QList<QString > vars, modelica_integer numberOfIntervals, modelica_boolean removeDescription);
  QList<QString > readSimulationResultVars(QString fileName, modelica_boolean readParameters, modelica_boolean openmodelicaStyle);
  modelica_integer readSimulationResultSize(QString fileName);
  modelica_boolean plotAll(modelica_boolean externalWindow, QString fileName, QString title, QString grid, modelica_boolean logX, modelica_boolean logY, QString xLabel, QString yLabel, QList<modelica_real > xRange, QList<modelica_real > yRange, modelica_real curveWidth, modelica_integer curveStyle, QString legendPosition, QString footer, modelica_boolean autoScale, modelica_boolean forceOMPlot);
  QList<QString > getAllSubtypeOf(QString parentClass, QString class_, modelica_boolean qualified, modelica_boolean includePartial, modelica_boolean sort);
  QList<QString > getPackages(QString class_);
  QList<QString > getUsedClassNames(QString className);
  QList<QString > getClassNames(QString class_, modelica_boolean recursive, modelica_boolean qualified, modelica_boolean sort, modelica_boolean builtin, modelica_boolean showProtected, modelica_boolean includeConstants);
  modelica_boolean setClassComment(QString class_, QString filename);
  modelica_boolean isShortDefinition(QString class_);
  modelica_boolean setSourceFile(QString class_, QString filename);
  QString getSourceFile(QString class_);
  modelica_boolean copyClass(QString className, QString newClassName, QString withIn);
  modelica_boolean moveClassToBottom(QString className);
  modelica_boolean moveClassToTop(QString className);
  modelica_boolean moveClass(QString className, modelica_integer offset);
  QList<QString > reduceTerms(QString className, modelica_real startTime, modelica_real stopTime, modelica_integer numberOfIntervals, modelica_real tolerance, QString method, QString fileNamePrefix, QString options, QString outputFormat, QString variableFilter, QString cflags, QString simflags, QString labelstoCancel);
  QList<QString > buildLabel(QString className, modelica_real startTime, modelica_real stopTime, modelica_integer numberOfIntervals, modelica_real tolerance, QString method, QString fileNamePrefix, QString options, QString outputFormat, QString variableFilter, QString cflags, QString simflags);
  modelica_boolean buildEncryptedPackage(QString className, modelica_boolean encrypt);
  QString buildModelFMU(QString className, QString version, QString fmuType, QString fileNamePrefix, QList<QString > platforms, modelica_boolean includeResources);
  QString translateModelFMU(QString className, QString version, QString fmuType, QString fileNamePrefix, modelica_boolean includeResources);
  QString importFMUModelDescription(QString filename, QString workdir, modelica_integer loglevel, modelica_boolean fullPath, modelica_boolean debugLogging, modelica_boolean generateInputConnectors, modelica_boolean generateOutputConnectors);
  QString importFMU(QString filename, QString workdir, modelica_integer loglevel, modelica_boolean fullPath, modelica_boolean debugLogging, modelica_boolean generateInputConnectors, modelica_boolean generateOutputConnectors);
  QList<QList<QString > > getLoadedLibraries();
  QString uriToFilename(QString uri);
  modelica_boolean rewriteBlockCall(QString className, QString inDefs);
  modelica_boolean generateVerificationScenarios(QString path);
  modelica_boolean inferBindings(QString path);
  modelica_boolean exportToFigaro(QString path, QString directory, QString database, QString mode, QString options, QString processor);
  QString listFile(QString class_, modelica_boolean nestedClasses);
  QString stringReplace(QString str, QString source, QString target);
  QList<QString > stringSplit(QString string, QString token);
  QList<QString > strtok(QString string, QString token);
  QList<QString > listVariables();
  QList<QString > getDerivedUnits(QString baseUnit);
  typedef struct {
    modelica_boolean unitsCompatible;
    modelica_real scaleFactor;
    modelica_real offset;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(unitsCompatible ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(QString::number(scaleFactor));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(offset));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } convertUnits_res;
  convertUnits_res convertUnits(QString s1, QString s2);
  typedef struct {
    modelica_boolean success;
    QString xmlfileName;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(success ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append("\"" + xmlfileName + "\"");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } dumpXMLDAE_res;
  dumpXMLDAE_res dumpXMLDAE(QString className, QString translationLevel, modelica_boolean addOriginalIncidenceMatrix, modelica_boolean addSolvingInfo, modelica_boolean addMathMLCode, modelica_boolean dumpResiduals, QString fileNamePrefix, QString rewriteRulesFile);
  QString translateGraphics(QString className);
  modelica_boolean save(QString className);
  modelica_boolean saveTotalModel(QString fileName, QString className, modelica_boolean stripAnnotations, modelica_boolean stripComments);
  modelica_boolean saveModel(QString fileName, QString className);
  modelica_boolean deleteFile(QString fileName);
  modelica_boolean loadModel(QString className, QList<QString > priorityVersion, modelica_boolean notify, QString languageStandard, modelica_boolean requireExactVersion);
  modelica_boolean generateCode(QString className);
  QString runOpenTURNSPythonScript(QString pythonScriptFile);
  QString buildOpenTURNSInterface(QString className, QString pythonTemplateFile, modelica_boolean showFlatModelica);
  QString instantiateModel(QString className);
  QString checkAllModelsRecursive(QString className, modelica_boolean checkProtected);
  QString checkModel(QString className);
  modelica_boolean remove(QString path);
  modelica_boolean copy(QString source, QString destination);
  modelica_boolean mkdir(QString newDirectory);
  QString cd(QString newWorkingDirectory);
  QString getAstAsCorbaString(QString fileName);
  QString getLanguageStandard();
  modelica_boolean getOrderConnections();
  modelica_boolean getShowAnnotations();
  modelica_boolean setShowAnnotations(modelica_boolean show);
  modelica_integer getDefaultOpenCLDevice();
  modelica_integer getVectorizationLimit();
  modelica_boolean setNoSimplify(modelica_boolean noSimplify);
  modelica_boolean getNoSimplify();
  QString getAnnotationVersion();
  QString getClassesInModelicaPath();
  modelica_boolean echo(modelica_boolean setEcho);
  QString runScript(QString fileName);
  modelica_boolean clearMessages();
  typedef struct {
    modelica_integer numMessages;
    modelica_integer numErrors;
    modelica_integer numWarnings;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(QString::number(numMessages));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(numErrors));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(numWarnings));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } countMessages_res;
  countMessages_res countMessages();
  QString getMessagesString();
  QString getErrorString(modelica_boolean warningsAsErrors);
  QString readFileNoNumeric(QString fileName);
  modelica_integer alarm(modelica_integer seconds);
  modelica_boolean compareFiles(QString file1, QString file2);
  modelica_boolean compareFilesAndMove(QString newFile, QString oldFile);
  modelica_boolean writeFile(QString fileName, QString data, modelica_boolean append);
  QString readFile(QString fileName);
  typedef struct {
    modelica_boolean success;
    modelica_real fileSize;
    modelica_real mtime;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append(success ? "true" : "false");
      resultBuffer.append(",");
      resultBuffer.append(QString::number(fileSize));
      resultBuffer.append(",");
      resultBuffer.append(QString::number(mtime));
      resultBuffer.append(")");
      return resultBuffer;
    }
  } stat_res;
  stat_res stat(QString fileName);
  QString getVersion(QString cl);
  modelica_boolean clearCommandLineOptions();
  typedef struct {
    QList<QString > validOptions;
    QString mainDescription;
    QList<QString > descriptions;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append("{");
      int validOptions_i = 0;
      foreach(QString validOptions_elt, validOptions) {
        if (validOptions_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + validOptions_elt + "\"");
        validOptions_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(",");
      resultBuffer.append("\"" + mainDescription + "\"");
      resultBuffer.append(",");
      resultBuffer.append("{");
      int descriptions_i = 0;
      foreach(QString descriptions_elt, descriptions) {
        if (descriptions_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + descriptions_elt + "\"");
        descriptions_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getConfigFlagValidOptions_res;
  getConfigFlagValidOptions_res getConfigFlagValidOptions(QString flag);
  QList<QString > getCommandLineOptions();
  modelica_boolean setCommandLineOptions(QString option);
  typedef struct {
    QList<QString > allChoices;
    QList<QString > allComments;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append("{");
      int allChoices_i = 0;
      foreach(QString allChoices_elt, allChoices) {
        if (allChoices_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allChoices_elt + "\"");
        allChoices_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(",");
      resultBuffer.append("{");
      int allComments_i = 0;
      foreach(QString allComments_elt, allComments) {
        if (allComments_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allComments_elt + "\"");
        allComments_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getAvailableTearingMethods_res;
  getAvailableTearingMethods_res getAvailableTearingMethods();
  QString getTearingMethod();
  typedef struct {
    QList<QString > allChoices;
    QList<QString > allComments;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append("{");
      int allChoices_i = 0;
      foreach(QString allChoices_elt, allChoices) {
        if (allChoices_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allChoices_elt + "\"");
        allChoices_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(",");
      resultBuffer.append("{");
      int allComments_i = 0;
      foreach(QString allComments_elt, allComments) {
        if (allComments_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allComments_elt + "\"");
        allComments_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getAvailableIndexReductionMethods_res;
  getAvailableIndexReductionMethods_res getAvailableIndexReductionMethods();
  QString getIndexReductionMethod();
  typedef struct {
    QList<QString > allChoices;
    QList<QString > allComments;
    QString toString() {
      QString resultBuffer = "(";
      resultBuffer.append("{");
      int allChoices_i = 0;
      foreach(QString allChoices_elt, allChoices) {
        if (allChoices_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allChoices_elt + "\"");
        allChoices_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(",");
      resultBuffer.append("{");
      int allComments_i = 0;
      foreach(QString allComments_elt, allComments) {
        if (allComments_i) {
          resultBuffer.append(",");
        }
        resultBuffer.append("\"" + allComments_elt + "\"");
        allComments_i++;
      }
      resultBuffer.append("}");
      resultBuffer.append(")");
      return resultBuffer;
    }
  } getAvailableMatchingAlgorithms_res;
  getAvailableMatchingAlgorithms_res getAvailableMatchingAlgorithms();
  QString getMatchingAlgorithm();
  modelica_boolean clearDebugFlags();
  modelica_boolean disableNewInstantiation();
  modelica_boolean enableNewInstantiation();
  modelica_boolean setCompilerFlags(QString compilerFlags);
  QString getModelicaPath();
  modelica_boolean setModelicaPath(QString modelicaPath);
  QString getInstallationDirectoryPath();
  modelica_boolean setInstallationDirectoryPath(QString installationDirectoryPath);
  modelica_boolean setEnvironmentVar(QString var, QString value);
  QString getEnvironmentVar(QString var);
  QString getTempDirectoryPath();
  modelica_boolean setTempDirectoryPath(QString tempDirectoryPath);
  modelica_boolean setPlotCommand(QString plotCommand);
  modelica_boolean setCompileCommand(QString compileCommand);
  QString getCompileCommand();
  modelica_boolean setCompilerPath(QString compilerPath);
  modelica_boolean verifyCompiler();
  modelica_boolean setCXXCompiler(QString compiler);
  QString getCXXCompiler();
  QString getCFlags();
  modelica_boolean setCFlags(QString inString);
  modelica_boolean setCompiler(QString compiler);
  QString getCompiler();
  modelica_boolean setLinkerFlags(QString linkerFlags);
  QString getLinkerFlags();
  modelica_boolean setLinker(QString linker);
  QString getLinker();
  modelica_boolean generateSeparateCodeDependenciesMakefile(QString filename, QString directory, QString suffix);
  QList<QString > generateSeparateCodeDependencies(QString stampSuffix);
  modelica_boolean generateSeparateCode(QString className, modelica_boolean cleanCache);
  modelica_boolean generateJuliaHeader(QString fileName);
  modelica_boolean generateHeader(QString fileName);
  modelica_boolean clearVariables();
  modelica_boolean clearProgram();
  modelica_boolean clear();
  QString help(QString topic);
  modelica_boolean saveAll(QString fileName);
  QList<modelica_integer > system_parallel(QList<QString > callStr, modelica_integer numThreads);
  modelica_integer system(QString callStr, QString outputFile);
  QList<QString > loadFileInteractive(QString filename, QString encoding);
  QList<QString > loadFileInteractiveQualified(QString filename, QString encoding);
  QList<QString > parseFile(QString filename, QString encoding);
  QList<QString > parseString(QString data, QString filename);
  modelica_boolean loadString(QString data, QString filename, QString encoding, modelica_boolean merge);
  modelica_boolean reloadClass(QString name, QString encoding);
  modelica_boolean loadEncryptedPackage(QString fileName, QString workdir);
  QList<QString > parseEncryptedPackage(QString fileName, QString workdir);
  modelica_boolean loadFiles(QList<QString > fileNames, QString encoding, modelica_integer numThreads);
  modelica_boolean loadFile(QString fileName, QString encoding, modelica_boolean uses);
signals:
  void logCommand(QString command, QTime *commandTime);
  void logResponse(QString command, QString response, QTime *responseTime);
  void throwException(QString exception);
};