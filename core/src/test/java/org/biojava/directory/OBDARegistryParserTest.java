/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.directory;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;

/**
 * <code>OBDARegistryParserTest</code> tests that registry config
 * files are correctly parsed.
 *
 * @author Keith James
 */
public class OBDARegistryParserTest extends TestCase
{
    public OBDARegistryParserTest(String name)
    {
        super(name);
    }

    public void testParseRegistry() throws Exception
    {
        RegistryConfiguration regConf = null;
        BufferedReader br = null;
        String locator = "<locator>";

        try
        {
            br = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("/org/biojava/directory/seqdatabase.ini")));
            regConf = OBDARegistryParser.parseRegistry(br, locator);
        }
        finally
        {
            if (br != null)
                br.close();
        }

        assertNotNull(regConf);
        assertEquals(locator, regConf.getConfigLocator());

        Map conf = regConf.getConfiguration();
        assertTrue(conf.containsKey("databank_0"));
        assertTrue(conf.containsKey("databank_1"));

        List dbConfigs = (List) conf.get("databank_0");
        assertEquals(1, dbConfigs.size());
 
        Map dbConfig = (Map) dbConfigs.get(0);
        assertEquals("<protocol-type 0>", dbConfig.get("protocol"));
        assertEquals("<location-string 0>", dbConfig.get("location"));

        dbConfigs = (List) conf.get("databank_1");
        assertEquals(2, dbConfigs.size());

        dbConfig = (Map) dbConfigs.get(0);
        assertEquals("<protocol-type 1a>", dbConfig.get("protocol"));
        assertEquals("<location-string 1a>", dbConfig.get("location"));

        dbConfig = (Map) dbConfigs.get(1);
        assertEquals("<protocol-type 1b>", dbConfig.get("protocol"));
        assertEquals("<location-string 1b>", dbConfig.get("location"));
    }
}
