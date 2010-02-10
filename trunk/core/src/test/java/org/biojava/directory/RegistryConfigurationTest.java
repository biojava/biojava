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

import java.util.HashMap;
import java.util.Map;

import junit.framework.TestCase;

/**
 * <code>RegistryConfigurationTest</code> tests basic behaviour of
 * <code>RegistryConfiguration</code>s. Note that the OBDA
 * implementation maps a database identifier string to a
 * <code>List</code> of <code>Map</code>s which acts as a sequence of
 * fallback mappings. The test maps the database identifier string
 * directly to config <code>Map</code>s for simplicity.
 *
 * @author Keith James
 */
public class RegistryConfigurationTest extends TestCase
{
    protected Map confParams0;
    protected Map confParams1;
    protected Map confParams2;
    protected String locator0;
    protected String locator1;
    protected String locator2;

    public RegistryConfigurationTest(String name)
    {
        super(name);
    }

    protected void setUp()
    {
        confParams0 = new HashMap();
        confParams0.put("protocol", "<protocol-type 0>");
        confParams0.put("location", "<location-string 0>");
        confParams0.put("parameter-x", "<value 0>");

        confParams1 = new HashMap();
        confParams1.put("protocol", "<protocol-type 1>");
        confParams1.put("location", "<location-string 1>");
        confParams1.put("parameter-y", "<value 1>");

        confParams2 = new HashMap();
        confParams2.put("protocol", "<protocol-type 2>");
        confParams2.put("location", "<location-string 2>");
        confParams2.put("parameter-x", "<value 2>");

        locator0 = "<locator0>";
        locator1 = "<locator1>";
        locator2 = "<locator2>";
    }

    public void testSimpleRegistry() throws Exception
    {
        Map conf0 = new HashMap();
        conf0.put("databank_0", confParams0);

        RegistryConfiguration simple =
            new RegistryConfiguration.Impl(locator0, conf0);

        assertEquals(locator0, simple.getConfigLocator());
        assertEquals(conf0, simple.getConfiguration());
    }

    public void testCompositeRegistry() throws Exception
    {
        RegistryConfiguration.Composite composite =
            new RegistryConfiguration.Composite();

        Map conf0 = new HashMap();
        conf0.put("databank_0", confParams0);
        RegistryConfiguration simple0 =
            new RegistryConfiguration.Impl(locator0, conf0);

        Map conf1 = new HashMap();
        conf1.put("databank_0", confParams1);
        RegistryConfiguration simple1 =
            new RegistryConfiguration.Impl(locator1, conf1);

        Map conf2 = new HashMap();
        conf2.put("databank_0", confParams2);
        RegistryConfiguration simple2 =
            new RegistryConfiguration.Impl(locator2, conf2);

        // Initial config
        composite.addTopConfig(simple0);
        assertEquals(locator0, composite.getConfigLocator());
        assertEquals("<protocol-type 0>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("protocol"));
        assertEquals("<location-string 0>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("location"));

        // This should be obscured by the initial config
        composite.addBottomConfig(simple1);
        assertEquals(locator0 + "::" + locator1,
                     composite.getConfigLocator());
        assertEquals("<protocol-type 0>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("protocol"));
        assertEquals("<location-string 0>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("location"));

        // This should obscure previous configs
        composite.addTopConfig(simple2);
        assertEquals(locator2 + "::" + locator0 + "::" + locator1,
                     composite.getConfigLocator());
        assertEquals("<protocol-type 2>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("protocol"));
        assertEquals("<location-string 2>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("location"));

        // Parameter present in top and middle; should get top value
        assertEquals("<value 2>",
                     ((Map) composite.getConfiguration().get("databank_0")).get("parameter-x"));

        // Parameter present in bottom only; should get null
        assertNull(((Map) composite.getConfiguration().get("databank_0")).get("parameter-y"));
    }
}
