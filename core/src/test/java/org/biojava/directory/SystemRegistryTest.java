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

import java.util.List;

import junit.framework.TestCase;

/**
 * <code>SystemRegistryTest</code> tests that a registry can be
 * created and that it will look in the correct places for a
 * configuration file.
 *
 * @author Keith James
 */
public class SystemRegistryTest extends TestCase
{
    public SystemRegistryTest(String name)
    {
        super(name);
    }

    public void testInstance() throws Exception
    {
        Registry registry = SystemRegistry.instance();
        assertNotNull(registry);
    }

    public void testGetRegistryPath()
    {
        List registryPath = SystemRegistry.getRegistryPath();
        int numElements = registryPath.size();

        if (numElements < 2 || numElements > 3)
            fail("Expected 2 or 3 registry path elements but found "
                 + numElements);

        int offset = 0;

        if (numElements == 3)
        {
            assertEquals("file://"
                         + System.getProperty("user.home")
                         + "/.bioinformatics/seqdatabase.ini",
                         registryPath.get(0));
            offset++;
        }

        assertEquals("file:///etc/bioinformatics/seqdatabase.ini",
                     registryPath.get(0 + offset));
        assertEquals("http://www.open-bio.org/registry/seqdatabase.ini",
                     registryPath.get(1 + offset));
    }
}
