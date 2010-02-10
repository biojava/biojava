package org.biojava.bio.symbol;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;

/**
 * Tests that serilization works as advertised.
 *
 * @author Thomas Down
 * @author Mark Schreiber
 * @since 1.3
 */

public class SymbolSerializationTest extends TestCase {
    public SymbolSerializationTest(String name){
        super(name);
    }

    private void doSymbolTest(Symbol s)
        throws Exception
    {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream((os));
        oos.writeObject(s);
        oos.flush();
        oos.close();

        ObjectInputStream ois = new ObjectInputStream(
                new ByteArrayInputStream(os.toByteArray()));
        Symbol s2 = (Symbol) ois.readObject();
        ois.close();

        assertTrue(s == s2);
    }

    public void testAtomic()
        throws Exception
    {
        doSymbolTest(DNATools.t());
    }

    public void testAmbiguous()
        throws Exception
    {
        doSymbolTest(DNATools.n());
    }
    
    
    public void testSpecialGap() throws Exception{
        doSymbolTest(AlphabetManager.getGapSymbol());
    }
    
    public void testCompoundGap() throws Exception{
        List alphas = Arrays.asList(new Alphabet[]{DNATools.getDNA(), DNATools.getDNA()});
        doSymbolTest(AlphabetManager.getGapSymbol(alphas));
    }
    
    public void testNormalGap() throws Exception{
        doSymbolTest(
                DNATools.getDNA().getGapSymbol());
    }
}
