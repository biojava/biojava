package org.biojava3.core.sequence.location;

import java.util.Arrays;
import java.util.Collection;
import org.biojava3.core.sequence.DataSource;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.CompoundSet;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

/**
 *
 * @author Jacek Grzebyta
 */
@RunWith(Parameterized.class)
public class TargetedLocationParserTest {
    
    private Data request;
    
    public static class Data {
        private String gi;
        private String Insdc;
        private CompoundSet<?> compound;
        private CompoundSet<?> originType;

        /**
         * Parser data input. Based on that input it should be able to identity the origin and wanted target
         * @param gi origin GI number
         * @param originType origin compound set type
         * @param Insdc string with INSDC notation
         * @param compound wanted compound type {@see CompoundSet}
         */
        public Data(String gi, CompoundSet originType, String Insdc, CompoundSet compound) {
            this.gi = gi;
            this.originType = originType;
            this.Insdc = Insdc;
            this.compound = compound;
        }
    };
    
    @Parameterized.Parameters
    public static Collection<Data[]> getLocations() throws Exception {
        

        Data[][] out = new Data[][]{
            {new Data("7525057", AminoAcidCompoundSet.getAminoAcidCompoundSet(), 
                    "join(complement(NC_000932.1:69611..69724),NC_000932.1:139856..140087,NC_000932.1:140625..140650)", DNACompoundSet.getDNACompoundSet())},
            
            {new Data("7525059", AminoAcidCompoundSet.getAminoAcidCompoundSet(), 
                    "NC_000932.1:72371..73897", DNACompoundSet.getDNACompoundSet())},
            
            {new Data("7525073", DNACompoundSet.getDNACompoundSet() ,
                    "complement(NC_000932.1:84005..84283)", DNACompoundSet.getDNACompoundSet())},
            
            
            {new Data("7525012", DNACompoundSet.getDNACompoundSet(), 
                    "complement(9938..11461)", DNACompoundSet.getDNACompoundSet())}
        };
        
        return Arrays.asList(out);
    }

    public TargetedLocationParserTest(Data request) {
        this.request = request;
    }
    
    
    @Test
    public void locationTest() throws Exception {
        
        InsdcParser parser = new InsdcParser(DataSource.GENBANK);
        Location loc = parser.parse(request.Insdc);
        
        Assert.assertNotNull(loc);
        if (loc.isComplex()) {
            Assert.assertFalse(loc.getSubLocations().isEmpty());
        } else {
        }
    }
}
