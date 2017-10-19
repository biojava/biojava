package org.biojava.nbio.core.sequence.io.embl;


import java.io.*;
import java.util.LinkedList;
import java.util.List;

public class EmblParser {

    private File file;
    private EmblId emblId = new EmblId();
    private EmblReference emblReference = new EmblReference();
    private List<String> accessionNumber = new LinkedList<>();
    private String ProjectIdentifier;
    private String OrGanelle;
    private String createdDate;
    private String lastUpdatedDate;
    private String sequenceDescription;
    private List<String> keyword = new LinkedList<>();
    private String organismSpecies;
    private String organismClassification;
    private String databaseCrossReference;
    private String assemblyHeader;
    private String assemblyInformation;
    private String CON;
    private String sequenceHeader;
    private StringBuilder sequence;

    public EmblId getEmblId() {
        return emblId;
    }

    public EmblReference getEmblReference() {
        return emblReference;
    }

    public List<String> getAccessionNumber() {
        return accessionNumber;
    }

    public String getProjectIdentifier() {
        return ProjectIdentifier;
    }

    public String getCreatedDate() {
        return createdDate;
    }

    public String getLastUpdatedDate() {
        return lastUpdatedDate;
    }

    public String getSequenceDescription() {
        return sequenceDescription;
    }

    public List<String> getKeyword() {
        return keyword;
    }

    public String getOrganismSpecies() {
        return organismSpecies;
    }

    public String getOrganismClassification() {
        return organismClassification;
    }


    public String getDatabaseCrossReference() {
        return databaseCrossReference;
    }

    public String getAssemblyHeader() {
        return assemblyHeader;
    }


    public String getAssemblyInformation() {
        return assemblyInformation;
    }


    public String getCON() {
        return CON;
    }


    public String getSequenceHeader() {
        return sequenceHeader;
    }


    public StringBuilder getSequence() {
        return sequence;
    }


    public EmblParser(File file) {
        setFile(file);
    }

    public void setFile(File file) {
        if (file == null)
            throw new NullPointerException("file can't be null");
        this.file = file;
    }

    public void parse() {
        try (FileReader fileReader = new FileReader(file)) {
            String line = "";
            String lineIdentifier;
            String lineInfo;
            try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
                while (bufferedReader.readLine() != null)
                    line = bufferedReader.readLine();
                lineInfo = line.substring(2);
                lineIdentifier = line.substring(0, 2);
                if (lineIdentifier.equals("ID"))
                    populateID(line);
                else if (lineIdentifier.equals("AC"))
                    populateAccessionNumber(line);
                else if (lineIdentifier.equals("DT") && line.contains("Created"))
                    createdDate = lineInfo;
                else if (lineIdentifier.equals("DT") && line.contains("updated"))
                    lastUpdatedDate = lineInfo;
                else if (lineIdentifier.equals("DE"))
                    sequenceDescription = lineInfo;
                else if (lineIdentifier.equals("KW"))
                    keyword.add(lineInfo);
                else if (lineIdentifier.equals("OS"))
                    organismSpecies = lineInfo;
                else if (lineIdentifier.equals("OC"))
                    organismClassification = lineInfo;
                else if (lineIdentifier.equals("OG"))
                    OrGanelle = lineInfo;
                else if (lineIdentifier.equals("RN") || lineIdentifier.equals("RP") || lineIdentifier.equals("RX")
                        || lineIdentifier.equals("RG") || lineIdentifier.equals("RA")
                        || lineIdentifier.equals("RT") || lineIdentifier.equals("RL"))
                    emblReferencePopulating(lineIdentifier, lineInfo);


            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void emblReferencePopulating(String lineIdentifier, String lineInfo) {
        if (lineIdentifier.equals("RN"))
            emblReference.setReferenceNumber(lineInfo);
        else if (lineIdentifier.equals("RP"))
            emblReference.setReferencePosition(lineInfo);
        else if (lineIdentifier.equals("RX"))
            emblReference.setReferenceCrossReference(lineInfo);
        else if (lineIdentifier.equals("RG"))
            emblReference.setReferenceGroup(lineInfo);
        else if (lineIdentifier.equals("RA"))
            emblReference.setReferenceAuthor(lineInfo);
        else if (lineIdentifier.equals("RT"))
            emblReference.setReferenceTitle(lineInfo);
        else if (lineIdentifier.equals("RL"))
            emblReference.setReferenceLocation(lineInfo);
    }

    private void populateAccessionNumber(String line) {
        accessionNumber.add(line);
    }

    private void populateID(String line) {
        line.replace(",", "");
        String[] strings = line.split(" ");
        emblId.setPrimaryAccession(strings[1]);
        emblId.setSequenceVersion(strings[2]);
        emblId.setTopology(strings[3]);
        emblId.setMoleculeType(strings[4]);
        emblId.setDataClass(strings[5]);
        emblId.setTaxonomicDivision(strings[6]);
        emblId.setSequenceLength(strings[7]);
    }


}
