package org.biojava.nbio.core.sequence;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TaxonomyIDTest {

    @Test
    void createTaxonomyID(){
        TaxonomyID tId = new TaxonomyID("abc1", DataSource.GENBANK);
        assertEquals("abc1", tId.getID());
        assertEquals(DataSource.GENBANK, tId.getDataSource());
    }

}