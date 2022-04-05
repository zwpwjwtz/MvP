#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 20:10:13 2022

Python script for parsing spectrum files of mzXML format
Required modules: base64, hashlib, re, struct, xml
"""


# Utility functions
def requireModule(name: str, alias = ''):
    if type(alias) != str or len(alias) == 0:
        alias = name
    try:
        globals()[alias] = __import__(name, fromlist = [None])
    except:
        print('Failed importing module "' + name + '" that is necessary \
for handling spectrum data. Please make sure you have installed it.')
        __import__('sys').exit(1)
    return True

def decodeBase64PeakList(peakString, precision):
    requireModule('base64')
    requireModule('struct')
    
    peakString = base64.b64decode(peakString)
    width = int(precision / 8)
    peakCount = int(len(peakString) / width)
    return [list(map(lambda i: 
                        struct.unpack('!f', peakString[i * 4 : i * 4 + 4])[0],
                     range(0, peakCount - 1, 2))),
            list(map(lambda i: 
                        struct.unpack('!f', peakString[i * 4 + 4 : i * 4 + 8])[0],
                     range(0, peakCount - 1, 2)))]

def encodeBase64PeakList(spectrum):
    requireModule('base64')
    requireModule('struct')
    
    peakList = list()
    for mz, intensity in zip(spectrum[0], spectrum[1]):
        peakList.append(struct.pack('!ff', mz, intensity))
    return str(base64.b64encode(b''.join(peak for peak in peakList)),
               encoding = 'utf-8')


# Handlers for the "TXT spectrum"
# This is the most commonly accept one, 
# with one m/z and its intensity per line seperated by 
# a separator (usually a whitespace)
# Using a comma (,) as separator will make it like a CSV file
def readTxtSpectrum(filename: str, fieldSeparator = ' ') -> list:
    mz = []
    intensities = []
    
    file = open(filename, 'rb')
    while True:
        line = str(file.readline(), encoding = 'utf-8')
        if line == '':
            break
        pos = line.find(fieldSeparator)
        if pos < 0:
            continue
        mz.append(float(line[0:pos]))
        intensities.append(float(line[pos + 1:(len(line) - 1)]))
    file.close()
    
    if len(mz) > 0:
        return [mz, intensities]
    else:
        return []
    
def writeTxtSpectrum(filename: str, 
                     spectrum: list, 
                     fieldSeparator = ' '):
    file = open(filename, 'wb')
    for mz, intensity in zip(spectrum[0], spectrum[1]):
        file.write('{}{}{}\n'.format(mz, fieldSeparator, intensity)
                             .encode('utf-8'))
    file.close()

# Handlers for the "mzXML spectrum"
def generateMZXMLHeader():
    return '<?xml version="1.0" encoding="ISO-8859-1"?>\n\
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2" \
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" \
xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2 \
http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_idx_3.2.xsd">\n\
<msRun scanCount="1" startTime="PT0S" endTime="PT0S">\n\
<dataProcessing>\n\
  <software type="conversion" name="unknown software name" version="2.16.2"/>\
  <processingOperation name="Conversion to mzXML"/>\n\
</dataProcessing>\n'

def readMZXML(filename: str) -> list:
    requireModule('re')
    requireModule('xml.etree.ElementTree', 'ElementTree')
    
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    namespace = re.match('\{.*\}', root.tag)
    if namespace.endpos > 0:
        namespace = namespace.group(0)
    else:
        namespace = ''
    scanList = list(root.iter('{}scan'.format(namespace)))
    peakList = list(map(lambda X: X.find('{}peaks'.format(namespace)), 
                        scanList))
    return list(map(decodeBase64PeakList, 
                    map(lambda X: X.text, peakList),
                    map(lambda X: int(X.attrib['precision']), peakList)))

def writeMZXML(filename: str, spectrumList: list):
    requireModule('hashlib')
    
    file = open(filename, 'wb')
    header = generateMZXMLHeader()
    body = '<scan num="1" scanType="Full" centroided="{}" msLevel="{}" \
    peaksCount="{}" polarity="{}" retentionTime="PT0S" \
lowMz="{}" highMz="{}" basePeakMz="{}" \
basePeakIntensity="{}" totIonCurrent="{}">\n\
<peaks precision="32" byteOrder="network" pairOrder="m/z-int">{}</peaks>\n\
</scan>\n</msRun>\n\
<index name="scan"><offset id="1">{}</offset></index>\n'
    tail1 = '<indexOffset>{}</indexOffset>\n'
    tail2 = '<sha1>{}</sha1>\n</mzXML>'
    
    spectrum = spectrumList[0]
    basePeakIndex = spectrum[1].index(max(spectrum[1]))
    body = body.format(0, # Centroided
                       1, # msLevel
                       len(spectrum[0]), # peakCount
                       '+', # polarity
                       min(spectrum[0]), # lowMz
                       max(spectrum[0]), # highMz
                       spectrum[0][basePeakIndex], # basePeakMz
                       spectrum[1][basePeakIndex], # basePeakIntensity
                       sum(spectrum[1]), # toIonCurrent
                       encodeBase64PeakList(spectrum), # peaks
                       len(header) # "offset" of the first scan
                      )
    tail1 = tail1.format(len(header) + len(body))
    tail2 = tail2.format(
            hashlib.sha1((header + body + tail1).encode('utf-8')).hexdigest())
    
    file.write((header + body + tail1 + tail2).encode('utf-8'))
    file.close()

# Read one spectrum from a data file by guessing its format
def readSpectrum(filename: str) -> list:
    spectrumList = []
    
	# Guess file type by suffix
    if filename.lower().endswith('.txt'):
	    spectrumList = readTxtSpectrum(filename)
    elif filename.lower().endswith('.csv'):
	    spectrumList = readTxtSpectrum(filename, fieldSeparator = ',')
    elif filename.lower().endswith('.mzxml'):
	    spectrumList = readMZXML(filename)
    else:
	    spectrumList = readTxtSpectrum(filename)
	
    if len(spectrumList) > 0 and type(spectrumList[0][0]) == list:
        # The file contains multiple spectra; return the fist spectrum only
        return spectrumList[0]
    else:
        return spectrumList

# Read a spectra file by guessing its format
def readSpectra(filename: str) -> list:
    # Guess file type by suffix
    if filename.lower().endswith('.txt'):
        return [readTxtSpectrum(filename)]
    elif filename.lower().endswith('.csv'):
        return [readTxtSpectrum(filename, fieldSeparator = ',')]
    elif filename.lower().endswith('.mzxml'):
        return readMZXML(filename)
    else:
        return []
