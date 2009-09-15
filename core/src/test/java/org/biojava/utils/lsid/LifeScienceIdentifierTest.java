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

package org.biojava.utils.lsid;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

/**
 * Test LifeScienceIdentifier.
 *
 * @author Michael Heuer
 */
public class LifeScienceIdentifierTest
    extends TestCase
{

    public LifeScienceIdentifierTest(String name)
    {
        super(name);
    }

    public void testValueOfParameters()
    {
        LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("authorityId",
                                                                   "namespaceId",
                                                                   "objectId",
                                                                   "revisionId");
        assertNotNull(lsid);
        assertTrue("authorityId".equals(lsid.getAuthorityId()));
        assertTrue("namespaceId".equals(lsid.getNamespaceId()));
        assertTrue("objectId".equals(lsid.getObjectId()));
        assertTrue("revisionId".equals(lsid.getRevisionId()));
    }
    public void testValueOfParametersNullRevisionId()
    {
        LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("authorityId",
                                                                   "namespaceId",
                                                                   "objectId",
                                                                   null);
        assertNotNull(lsid);
        assertTrue("authorityId".equals(lsid.getAuthorityId()));
        assertTrue("namespaceId".equals(lsid.getNamespaceId()));
        assertTrue("objectId".equals(lsid.getObjectId()));
        assertTrue(lsid.getRevisionId() == null);
    }
    public void testValueOfParametersEmptyRevisionId()
    {
        LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("authorityId",
                                                                   "namespaceId",
                                                                   "objectId",
                                                                   "");
        assertNotNull(lsid);
        assertTrue("authorityId".equals(lsid.getAuthorityId()));
        assertTrue("namespaceId".equals(lsid.getNamespaceId()));
        assertTrue("objectId".equals(lsid.getObjectId()));
        assertTrue("".equals(lsid.getRevisionId()));
    }
    public void testValueOfString()
    {
        try
        {
            LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("urn:lsid:authorityId:namespaceId:objectId:revisionId");
            
            assertNotNull(lsid);
            assertTrue("authorityId".equals(lsid.getAuthorityId()));
            assertTrue("namespaceId".equals(lsid.getNamespaceId()));
            assertTrue("objectId".equals(lsid.getObjectId()));
            assertTrue("revisionId".equals(lsid.getRevisionId()));
        }
        catch (LifeScienceIdentifierParseException e)
        {
            fail(e.getMessage());
        }
    }
    public void testValueOfStringNullRevisionId()
    {
        try
        {
            LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("urn:lsid:authorityId:namespaceId:objectId");
            
            assertNotNull(lsid);
            assertTrue("authorityId".equals(lsid.getAuthorityId()));
            assertTrue("namespaceId".equals(lsid.getNamespaceId()));
            assertTrue("objectId".equals(lsid.getObjectId()));
            assertTrue(lsid.getRevisionId() == null);
        }
        catch (LifeScienceIdentifierParseException e)
        {
            fail(e.getMessage());
        }
    }
    public void testValueOfStringEmptyRevisionId()
    {
        try
        {
            LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf("urn:lsid:authorityId:namespaceId:objectId:");
            
            assertNotNull(lsid);
            assertTrue("authorityId".equals(lsid.getAuthorityId()));
            assertTrue("namespaceId".equals(lsid.getNamespaceId()));
            assertTrue("objectId".equals(lsid.getObjectId()));
            assertTrue("".equals(lsid.getRevisionId()));
        }
        catch (LifeScienceIdentifierParseException e)
        {
            fail(e.getMessage());
        }
    }
    
    private boolean throwsParseException(String value)
    {
        try
        {
            LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf(value);
            lsid=lsid==null?null:lsid;//trick
        }
        catch (LifeScienceIdentifierParseException e)
        {
            return true;
        }
        
        return false;
    }
    
    public void testParseErrors()
    {
        List parseErrors = Arrays.asList(new String[] {"",
                                                       "urn:",
                                                       "urn:lsid",
                                                       "urn:lsid:",
                                                       "urn:lsid:authorityId",
                                                       "urn:lsid:authorityId:namespaceId",
                                                       "urn:lsid:authorityId:namespaceId:objectId:revisionId:extraAttribute",
                                                       "authority:",
                                                       "authority:namespace:",
                                                       "authority:namespace:objectId:"
        });
        
        for (Iterator i = parseErrors.iterator(); i.hasNext(); )
        {
            String value = (String) i.next();
            if (! (throwsParseException(value)))
            {
                fail("expected LifeScienceIdentifierParseException for: " + value);
            }
        }
    }

    private boolean throwsNullPointerException(String authorityId,
                                               String namespaceId,
                                               String objectId)
    {
        try
        {
            LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf(authorityId,
                                                                       namespaceId,
                                                                       objectId);
            lsid=lsid==null?null:lsid;//trick
        }
        catch (NullPointerException e)
        {
            return true;
        }
        
        return false;
    }

    public void testValueOfThrowsNullPointerException()
    {
        if (! (throwsNullPointerException(null, "namespaceId", "objectId")))
            fail("expected NullPointerException for null, namespaceId, objectId");

        if (! (throwsNullPointerException("authorityId", null, "objectId")))
            fail("expected NullPointerException for authorityId, null, objectId");

        if (! (throwsNullPointerException("authorityId", "namespaceId", null)))
            fail("expected NullPointerException for authorityId, namespaceId, null");
    }

    public void testObjectContract()
    {
        LifeScienceIdentifier lsid0 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "revisionId");
        
        LifeScienceIdentifier lsid1 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "revisionId");
        
        LifeScienceIdentifier lsid2 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "revisionId_1");
        
        LifeScienceIdentifier lsid3 = LifeScienceIdentifier.valueOf("authorityId_1",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "revisionId");
        
        LifeScienceIdentifier lsid4 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    null);
        
        LifeScienceIdentifier lsid5 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "");
        
        assertTrue(lsid0 != lsid1);
        assertTrue(lsid0.equals(lsid1));
        assertTrue(lsid0.hashCode() == lsid1.hashCode());
        assertTrue(lsid1 != lsid0);
        assertTrue(lsid1.equals(lsid0));
        assertTrue(lsid1.hashCode() == lsid0.hashCode());
        
        assertTrue(lsid0.equals(lsid2) == false);
        assertTrue(lsid0.hashCode() != lsid2.hashCode());
        assertTrue(lsid2.equals(lsid0) == false);
        assertTrue(lsid2.hashCode() != lsid0.hashCode());
        
        assertTrue(lsid0.equals(lsid3) == false);
        assertTrue(lsid0.hashCode() != lsid3.hashCode());
        assertTrue(lsid3.equals(lsid0) == false);
        assertTrue(lsid3.hashCode() != lsid0.hashCode());
        
        assertTrue(lsid0.equals(lsid4) == false);
        assertTrue(lsid0.hashCode() != lsid4.hashCode());
        assertTrue(lsid4.equals(lsid0) == false);
        assertTrue(lsid4.hashCode() != lsid0.hashCode());
        assertTrue(lsid0.equals(lsid5) == false);
        assertTrue(lsid0.hashCode() != lsid5.hashCode());
        assertTrue(lsid5.equals(lsid0) == false);
        assertTrue(lsid5.hashCode() != lsid0.hashCode());
        assertTrue(lsid4.equals(lsid5) == false);
        assertTrue(lsid4.hashCode() != lsid5.hashCode());
        assertTrue(lsid5.equals(lsid4) == false);
        assertTrue(lsid5.hashCode() != lsid4.hashCode());
        
        assertTrue("urn:lsid:authorityId:namespaceId:objectId:revisionId".equals(lsid0.toString()));
        assertTrue("urn:lsid:authorityId:namespaceId:objectId:revisionId".equals(lsid1.toString()));
        assertTrue("urn:lsid:authorityId:namespaceId:objectId:revisionId_1".equals(lsid2.toString()));
        assertTrue("urn:lsid:authorityId_1:namespaceId:objectId:revisionId".equals(lsid3.toString()));
        assertTrue("urn:lsid:authorityId:namespaceId:objectId".equals(lsid4.toString()));
        assertTrue("urn:lsid:authorityId:namespaceId:objectId:".equals(lsid5.toString()));
    }
    public void testSerialization()
    {
        LifeScienceIdentifier lsid0 = LifeScienceIdentifier.valueOf("authorityId",
                                                                    "namespaceId",
                                                                    "objectId",
                                                                    "revisionId");
        LifeScienceIdentifier lsid1 = null;
        
        ObjectOutputStream oos = null;
        ObjectInputStream ois = null;
        try 
        {
            //File f = File.createTempFile("lsid", ".ser");
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            oos = new ObjectOutputStream(new BufferedOutputStream(baos));
            oos.writeObject(lsid0);
            oos.flush();
            
            ois = new ObjectInputStream(new BufferedInputStream(new ByteArrayInputStream(baos.toByteArray())));
            lsid1 = (LifeScienceIdentifier) ois.readObject();
            
            assertTrue(lsid0 != lsid1);
            assertTrue(lsid0.equals(lsid1) == true);
        }
        catch (Exception e)
        {
            fail(e.getMessage());
        }
        finally
        {
            try
            {
                oos.close();
            }
            catch (Exception e1) {}
            try
            {
                ois.close();
            }
            catch (Exception e2) {}
        }
    }
}
