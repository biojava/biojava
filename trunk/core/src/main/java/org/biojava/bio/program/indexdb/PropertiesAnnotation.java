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

package org.biojava.bio.program.indexdb;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Map;
import java.util.Properties;

import org.biojava.bio.AbstractAnnotation;
import org.biojava.bio.BioError;
import org.biojava.utils.Commitable;

/**
 * <code>PropertiesAnnotation</code> creates an
 * <code>Annotation</code> view of a <code>Properties</code> bundle.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
class PropertiesAnnotation
    extends AbstractAnnotation implements Commitable {
    private Properties props;
    private File propsFile;

    /**
     * Creates a new <code>PropertiesAnnotation</code>.
     *
     * @param propsFile a <code>File</code> which should be a Java
     * <code>Properties</code> file.
     */
    public PropertiesAnnotation(File propsFile) {
        this.propsFile = propsFile;
        this.props = new Properties();

        if(propsFile.exists()) {
            try {
                props.load(new FileInputStream(propsFile));
            } catch (IOException ioe) {
                throw new BioError( "Assertion Failure: could not load properties",ioe);
            }
        }
    }

    public void commit() {
        try {
            props.store(new FileOutputStream(propsFile), "Meta-Data");
        } catch (IOException ioe) {
            try {
                rollback();
            } catch (BioError be) {
                throw new BioError("Catastrophic failure: could not roll back after failed commit",be);
            }
            throw new BioError("Could not commit");
        }
    }

    public void rollback() {
        if(propsFile.exists()) {
            try {
                props.load(new FileInputStream(propsFile));
            } catch (IOException ioe) {
                throw new BioError("Could not roll back");
            }
        } else {
            props.clear();
        }
    }

    protected Map getProperties() {
        return props;
    }

    protected boolean propertiesAllocated() {
        return true;
    }
}
