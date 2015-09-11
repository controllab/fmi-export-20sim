<?xml version="1.0"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- Mind the CDATA section: This is necessary to force output tags to use CDATA encapsulating -->
<xsl:output method="xml" version="1.0" encoding="ISO-8859-1" indent="yes"/>

<xsl:strip-space elements="*" />

<!-- External input given by the command line -->
<xsl:param name="SOURCEDIRECTORY"/> <!-- token that specifies the source file location -->

<!-- Make the number of variables globally available, and compute them (hopefully), once -->
<xsl:variable name="NUMBER_CONSTANTS" >
	<xsl:choose>
		<xsl:when test="/modelConfiguration/storage/tables/table[name='xx_C']">
			<xsl:value-of select="/modelConfiguration/storage/tables/table[name='xx_C']/count" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:variable name="NUMBER_PARAMETERS" >
	<xsl:choose>
		<xsl:when test="/modelConfiguration/storage/tables/table[name='xx_P']">
			<xsl:value-of select="/modelConfiguration/storage/tables/table[name='xx_P']/count" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:variable name="NUMBER_INITIAL_VALUES" >
	<xsl:choose>
		<xsl:when test="/modelConfiguration/storage/tables/table[name='xx_I']">
			<xsl:value-of select="/modelConfiguration/storage/tables/table[name='xx_I']/count" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:variable name="NUMBER_VARIABLES" >
	<xsl:choose>
		<xsl:when test="/modelConfiguration/storage/tables/table[name='xx_V']">
			<xsl:value-of select="/modelConfiguration/storage/tables/table[name='xx_V']/count" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:variable name="NUMBER_STATES" >
	<xsl:choose>
		<xsl:when test="/modelConfiguration/storage/tables/table[name='xx_s']">
			<xsl:value-of select="/modelConfiguration/storage/tables/table[name='xx_s']/count" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:variable name="NUMBER_RATES" >
	<xsl:value-of select="$NUMBER_STATES" />
</xsl:variable>

<!--  Template for root MCF element -->
<xsl:template match="/modelConfiguration">
	<xsl:element name="fmiModelDescription">
		<xsl:attribute name="fmiVersion">%FMIVERSION%</xsl:attribute>
		<xsl:attribute name="modelName">
			<xsl:value-of select="general/name" />
		</xsl:attribute>
%IF%%FMI1%
		<xsl:attribute name="modelIdentifier">
			<xsl:value-of select="general/name" />
		</xsl:attribute>
%ENDIF%
		<xsl:attribute name="guid">
			<xsl:text>{</xsl:text>
			<xsl:value-of select="document(concat($SOURCEDIRECTORY, '\GUID.xml'))/tokens/token[@name='GUID']" />
			<xsl:text>}</xsl:text>
		</xsl:attribute>
		<xsl:attribute name="generationTool">20-sim</xsl:attribute>
%IF%%FMI1%
		<xsl:attribute name="numberOfContinuousStates">
			<!-- Note: 20-sim doesn't distinguish between discrete and continuous states, so return the union of both sets -->
			<xsl:value-of select="$NUMBER_STATES" />
		</xsl:attribute>
%ENDIF%
		<xsl:attribute name="numberOfEventIndicators">0</xsl:attribute>

%IF%%FMI2%
		<xsl:attribute name="copyright">Controllab Products B.V.</xsl:attribute>
		<xsl:attribute name="license">-</xsl:attribute>
%ENDIF%

%IF%%FMI2%
		<xsl:element name="CoSimulation">
			<xsl:attribute name="modelIdentifier">
				<xsl:value-of select="general/name" />
			</xsl:attribute>
			<xsl:attribute name="needsExecutionTool">false</xsl:attribute>
			<xsl:attribute name="canHandleVariableCommunicationStepSize">false</xsl:attribute>
			<xsl:attribute name="canInterpolateInputs">true</xsl:attribute>
			<xsl:attribute name="maxOutputDerivativeOrder">0</xsl:attribute>
			<xsl:attribute name="canRunAsynchronuously">false</xsl:attribute>
			<xsl:attribute name="canBeInstantiatedOnlyOncePerProcess">true</xsl:attribute>
			<xsl:attribute name="canNotUseMemoryManagementFunctions">true</xsl:attribute>
			<xsl:attribute name="canGetAndSetFMUstate">false</xsl:attribute>
			<xsl:attribute name="canSerializeFMUstate">false</xsl:attribute>
			<xsl:attribute name="providesDirectionalDerivative">false</xsl:attribute>
		</xsl:element>
%ENDIF%

		<xsl:element name="DefaultExperiment">
			<xsl:attribute name="startTime"><xsl:value-of select="run/startTime" /></xsl:attribute>
			<xsl:attribute name="stopTime"><xsl:value-of select="run/finishTime" /></xsl:attribute>		
		</xsl:element>
		
		<xsl:element name="ModelVariables">
			<!-- Only real support currently -->
			<xsl:for-each select="modelVariables/modelVariable[string(type) = 'real']">
				<xsl:call-template name="modelVariable">
					<xsl:with-param name="modelvariable" select="."/>
					<xsl:with-param name="isArray" select="boolean(size and (number(size/rows) > 1 or number(size/columns) > 1)) "/>
					<xsl:with-param name="index" select="0"/>
				</xsl:call-template>
			</xsl:for-each>
		</xsl:element>
%IF%%FMI1%
		<xsl:element name="Implementation">
			<xsl:element name="CoSimulation_StandAlone">
				<xsl:element name="Capabilities">
					<xsl:attribute name="canHandleVariableCommunicationStepSize">false</xsl:attribute>
					<xsl:attribute name="canHandleEvents">false</xsl:attribute>
					<xsl:attribute name="canRejectSteps">false</xsl:attribute>
					<xsl:attribute name="canInterpolateInputs">true</xsl:attribute>
					<xsl:attribute name="maxOutputDerivativeOrder">0</xsl:attribute>
					<xsl:attribute name="canRunAsynchronuously">false</xsl:attribute>
					<xsl:attribute name="canSignalEvents">false</xsl:attribute>
					<xsl:attribute name="canBeInstantiatedOnlyOncePerProcess">true</xsl:attribute>
					<xsl:attribute name="canNotUseMemoryManagementFunctions">true</xsl:attribute>
				</xsl:element>
			</xsl:element>		
		</xsl:element>
%ENDIF%
%IF%%FMI2%
		<xsl:element name="ModelStructure">
			<xsl:if test="modelVariables/modelVariable[string(kind) = 'output']">
				<xsl:element name="Outputs">
					<xsl:for-each select="modelVariables/modelVariable[string(kind) = 'output']">
						<!-- Add 1-based index to the corresponding model variable -->
						<xsl:element name="Unknown">
							<xsl:attribute name="index">
								<xsl:value-of select="count(preceding-sibling::modelVariable) + 1" />
							</xsl:attribute>
						</xsl:element>	
					</xsl:for-each>				
				</xsl:element>			
			</xsl:if>			
		</xsl:element>
%ENDIF%
		
    </xsl:element>
</xsl:template>

<!-- Map each Model Configuration Model variable to a FMU ScalarVariable 
	@modelVariable	The model variable
	@isArray		Boolean indicating if the current variable is an array
	@index			The index within the array. Zero for scalar variables
-->
<xsl:template name="modelVariable">
	<xsl:param name="modelvariable"/>
	<xsl:param name="isArray"/>
	<xsl:param name="index"/>

	<xsl:variable name="variable_table" select="$modelvariable/storage/name" />
	<xsl:variable name="variable_offset" select="$modelvariable/storage/index" />

	<xsl:element name="ScalarVariable">
		<xsl:attribute name="name">
			<xsl:choose>
				<xsl:when test="$isArray"><xsl:value-of select="concat(translate($modelvariable/name,'\','.'), '[', $index, ']')" /></xsl:when>
				<xsl:otherwise><xsl:value-of select="translate($modelvariable/name,'\','.')" /></xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
		<xsl:attribute name="valueReference">
			<xsl:call-template name="valueReference">
				<xsl:with-param name="modelvariable" select="$modelvariable" />
				<xsl:with-param name="index" select="$index" />
			</xsl:call-template>
		</xsl:attribute>
		<xsl:if test="description and string-length(string($modelvariable/description)) > 0">
			<xsl:attribute name="description"><xsl:value-of select="$modelvariable/description" /></xsl:attribute>
		</xsl:if>
		<xsl:attribute name="variability">
			<xsl:choose>
%IF%%FMI1%
				<xsl:when test="string($modelvariable/kind)='parameter'">parameter</xsl:when>
%ENDIF%
%IF%%FMI2%
				<xsl:when test="string($modelvariable/kind)='parameter'">tunable</xsl:when>
				<xsl:when test="string($modelvariable/kind)='initial value'">fixed</xsl:when>
%ENDIF%
				<xsl:when test="string($modelvariable/kind)='constant'">constant</xsl:when>
				<xsl:when test="string(document(concat($SOURCEDIRECTORY, '\tokens.xml'))/tokens/token[@name='MODEL_IS_DISCRETE']) = 'XXTRUE'">discrete</xsl:when>
				<xsl:otherwise>continuous</xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
		<xsl:attribute name="causality">
			<xsl:choose>
%IF%%FMI2%
				<xsl:when test="string($modelvariable/kind)='parameter'">parameter</xsl:when>
				<xsl:when test="string($modelvariable/kind)='initial value'">parameter</xsl:when>
%ENDIF%
				<xsl:when test="string($modelvariable/kind)='input'">input</xsl:when>
				<xsl:when test="string($modelvariable/kind)='output'">output</xsl:when>
%IF%%FMI1%
				<xsl:otherwise>internal</xsl:otherwise>
%ENDIF%
%IF%%FMI2%
				<xsl:otherwise>local</xsl:otherwise>
%ENDIF%
			</xsl:choose>
		</xsl:attribute>
		<!--<xsl:if test="$modelvariable/preceding-sibling::modelVariable/storage[string(name) = $variable_table  and string(index) = $variable_offset] or $modelvariable/aliasOf">  Test if it's not the first usage of the underlying variable -->
%IF%%FMI1%
		<xsl:if test="$modelvariable/aliasOf">
			<xsl:attribute name="alias">alias</xsl:attribute>
		</xsl:if>		
%ENDIF%
		<!-- end of attributes assignment -->
	
		<xsl:element name="Real">
%IF%%FMI1%
			<xsl:if test="$modelvariable/value">
				<xsl:attribute name="start">
					<xsl:call-template name="GetArrayValue">
						<xsl:with-param name="stringarray" select="$modelvariable/value"/>
						<xsl:with-param name="index" select="$index"/>
					</xsl:call-template>
				</xsl:attribute>
				<!-- determine the 'fixed' attribute -->
				<xsl:attribute name="fixed">
					<xsl:choose>
						<xsl:when test="string($modelvariable/kind)='initial value'">true</xsl:when>
						<xsl:otherwise>false</xsl:otherwise>
					</xsl:choose>
				</xsl:attribute>
			</xsl:if>
%ENDIF%
%IF%%FMI2%
			<xsl:if test="string($modelvariable/kind)='parameter' or string($modelvariable/kind)='initial value' or string($modelvariable/kind)='input'">
				<xsl:attribute name="start">
					<xsl:call-template name="GetArrayValue">
						<xsl:with-param name="stringarray" select="$modelvariable/value"/>
						<xsl:with-param name="index" select="$index"/>
					</xsl:call-template>
				</xsl:attribute>
			</xsl:if>
%ENDIF%

		</xsl:element>
	</xsl:element>

	<!-- In case of an array, retrieve the next variable -->
	<xsl:if test="$isArray and number($index + 1) &lt; number($modelvariable/size/rows) * number($modelvariable/size/columns)">
		<xsl:call-template name="modelVariable">
			<xsl:with-param name="modelvariable" select="$modelvariable"/>
			<xsl:with-param name="isArray" select="$isArray"/>
			<xsl:with-param name="index" select="number($index + 1)"/>
		</xsl:call-template>
	</xsl:if>
</xsl:template>

<!-- Get the valueReference of a modelVariable -->
<xsl:template name="valueReference">
	<xsl:param name="modelvariable" />
	<xsl:param name="index" />

	<xsl:choose>
		<!-- Note that the data model (i.e. the storage array) defines the mapping to the FMU array, and that the table is not necessarily equal to the variable 'kind' -->
		<xsl:when test="string($modelvariable/storage/name)='xx_C'">
			<xsl:value-of select="$modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:when test="string($modelvariable/storage/name)='xx_P'">
			<xsl:value-of select="$NUMBER_CONSTANTS + $modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:when test="string($modelvariable/storage/name)='xx_I'">
			<xsl:value-of select="$NUMBER_CONSTANTS + $NUMBER_PARAMETERS + $modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:when test="string($modelvariable/storage/name)='xx_V'">
			<xsl:value-of select="$NUMBER_CONSTANTS + $NUMBER_PARAMETERS + $NUMBER_INITIAL_VALUES + $modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:when test="string($modelvariable/storage/name)='xx_s'">
			<xsl:value-of select="$NUMBER_CONSTANTS + $NUMBER_PARAMETERS + $NUMBER_INITIAL_VALUES + $NUMBER_VARIABLES + $modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:when test="string($modelvariable/storage/name)='xx_R'">
			<xsl:value-of select="$NUMBER_CONSTANTS + $NUMBER_PARAMETERS + $NUMBER_INITIAL_VALUES + $NUMBER_VARIABLES + $NUMBER_STATES + $modelvariable/storage/index + $index" />
		</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<!-- Get the value of a 20-sim stored array. 
	Recursive function that decreases the index until 0.
 
	@stringarray	Double colomn separate string, e.g. 0.0; 1.0; 1.20-sim
	@index			Index of the variable to be returned
	@return			Value at position index, or an empty object if nothing is found
-->
<xsl:template name="GetArrayValue">
	<xsl:param name="stringarray" />
	<xsl:param name="index" />
	
	<xsl:choose>
		<!-- if the ',' is not preceded by a ';' -->
		<xsl:when test="$index = 0 and contains($stringarray, ',') and string-length(substring-after($stringarray, ',')) &gt; string-length(substring-after($stringarray, ';'))">
			<xsl:value-of select="normalize-space(substring-before($stringarray, ','))" />
		</xsl:when>
		<xsl:when test="$index = 0 and contains($stringarray, ';')">
			<xsl:value-of select="normalize-space(substring-before($stringarray, ';'))" />
		</xsl:when>
		<xsl:when test="$index = 0">
			<xsl:value-of select="normalize-space($stringarray)" />
		</xsl:when>
		<!-- Matrix support. Make sure to return the row items if there are any before moving to the next column -->
		<xsl:when test="contains($stringarray, ',') and string-length(substring-after($stringarray, ',')) &gt; string-length(substring-after($stringarray, ';'))">
			<xsl:call-template name="GetArrayValue">
				<xsl:with-param name="stringarray" select="substring-after($stringarray, ',')"/>
				<xsl:with-param name="index" select="number($index - 1)"/>
			</xsl:call-template>
		</xsl:when>
		<xsl:when test="contains($stringarray, ';')">
			<xsl:call-template name="GetArrayValue">
				<xsl:with-param name="stringarray" select="substring-after($stringarray, ';')"/>
				<xsl:with-param name="index" select="number($index - 1)"/>
			</xsl:call-template>
		</xsl:when>
		<xsl:otherwise></xsl:otherwise>
	</xsl:choose>
</xsl:template>

</xsl:stylesheet> 

