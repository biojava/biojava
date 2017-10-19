package org.biojava.nbio.core.sequence.io.embl;


import java.io.*;

public class EmblParser {

    private File file;
    private EmblId id;
    private EmblReference emblReference;
    private String accessionNumber;
    private String ProjectIdentifier;
    private String createdDate;
    private String lastUpdatedDate;
    private String sequenceDescription;
    private String keyword;
    private String organismSpecies;
    private String organismClassification;
    private String databaseCrossReference;
    private String assemblyHeader;
    private String assemblyInformation;
    private String CON;
    private String sequenceHeader;
    private StringBuilder sequence;

    public EmblId getEmblId() {
        return id;
    }

    public EmblReference getEmblReference() {
        return emblReference;
    }


    public String getAccessionNumber() {
        return accessionNumber;
    }

    public void setAccessionNumber(String accessionNumber) {
        this.accessionNumber = accessionNumber;
    }

    public String getProjectIdentifier() {
        return ProjectIdentifier;
    }

    public void setProjectIdentifier(String projectIdentifier) {
        ProjectIdentifier = projectIdentifier;
    }

    public String getCreatedDate() {
        return createdDate;
    }

    public void setCreatedDate(String createdDate) {
        this.createdDate = createdDate;
    }

    public String getLastUpdatedDate() {
        return lastUpdatedDate;
    }

    public void setLastUpdatedDate(String lastUpdatedDate) {
        this.lastUpdatedDate = lastUpdatedDate;
    }

    public String getSequenceDescription() {
        return sequenceDescription;
    }

    public void setSequenceDescription(String sequenceDescription) {
        this.sequenceDescription = sequenceDescription;
    }

    public String getKeyword() {
        return keyword;
    }

    public void setKeyword(String keyword) {
        this.keyword = keyword;
    }

    public String getOrganismSpecies() {
        return organismSpecies;
    }

    public void setOrganismSpecies(String organismSpecies) {
        this.organismSpecies = organismSpecies;
    }

    public String getOrganismClassification() {
        return organismClassification;
    }

    public void setOrganismClassification(String organismClassification) {
        this.organismClassification = organismClassification;
    }

    public String getDatabaseCrossReference() {
        return databaseCrossReference;
    }

    public void setDatabaseCrossReference(String databaseCrossReference) {
        this.databaseCrossReference = databaseCrossReference;
    }

    public String getAssemblyHeader() {
        return assemblyHeader;
    }

    public void setAssemblyHeader(String assemblyHeader) {
        this.assemblyHeader = assemblyHeader;
    }

    public String getAssemblyInformation() {
        return assemblyInformation;
    }

    public void setAssemblyInformation(String assemblyInformation) {
        this.assemblyInformation = assemblyInformation;
    }

    public String getCON() {
        return CON;
    }

    public void setCON(String CON) {
        this.CON = CON;
    }

    public String getSequenceHeader() {
        return sequenceHeader;
    }

    public void setSequenceHeader(String sequenceHeader) {
        this.sequenceHeader = sequenceHeader;
    }

    public StringBuilder getSequence() {
        return sequence;
    }

    public void setSequence(StringBuilder sequence) {
        this.sequence = sequence;
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
            try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
                while (bufferedReader.readLine() != null)
                    line = bufferedReader.readLine();
                if (line.substring(0, 2).equals("ID"))
                    populateID(line);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void populateID(String line) {
        
    }


}
