// $Id: AnnotationParser.java,v 1.3 2009/11/03 19:16:34 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.io.parsers.phyloxml.phylogenydata;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Uri;

public class AnnotationParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new AnnotationParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private AnnotationParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        final Annotation annotation = new Annotation();
        if ( element.isHasAttribute( PhyloXmlMapping.ANNOTATION_TYPE_ATTR ) ) {
            annotation.setType( element.getAttribute( PhyloXmlMapping.ANNOTATION_TYPE_ATTR ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR ) ) {
            annotation.setEvidence( element.getAttribute( PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.ANNOTATION_REF_ATTR ) ) {
            annotation.setRef( element.getAttribute( PhyloXmlMapping.ANNOTATION_REF_ATTR ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.ANNOTATION_SOURCE_ATTR ) ) {
            annotation.setSource( element.getAttribute( PhyloXmlMapping.ANNOTATION_SOURCE_ATTR ) );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.ANNOTATION_DESC ) ) {
                annotation.setDesc( child_element.getValueAsString() );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.CONFIDENCE ) ) {
                annotation.setConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.URI ) ) {
                annotation.setUri( ( Uri ) UriParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.PROPERTY ) ) {
                if ( annotation.getProperties() == null ) {
                    annotation.setProperties( new PropertiesMap() );
                }
                annotation.getProperties()
                        .addProperty( ( Property ) PropertyParser.getInstance().parse( child_element ) );
            }
        }
        return annotation;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
