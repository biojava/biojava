// $Id: Annotation.java,v 1.24 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Annotation implements PhylogenyData {

    private String        _desc;
    private String        _type;
    private String        _source;
    private String        _ref;
    private String        _evidence;
    private Confidence    _confidence;
    private PropertiesMap _properties;
    private Uri           _uri;

    public Annotation() {
        init();
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getDesc() );
    }

    @Override
    public StringBuffer asText() {
        return new StringBuffer( getDesc() );
    }

    @Override
    public PhylogenyData copy() {
        final Annotation ann = new Annotation();
        if ( getConfidence() != null ) {
            ann.setConfidence( ( Confidence ) getConfidence().copy() );
        }
        else {
            ann.setConfidence( null );
        }
        ann.setType( new String( getType() ) );
        ann.setDesc( new String( getDesc() ) );
        ann.setEvidence( new String( getEvidence() ) );
        ann.setRef( new String( getRef() ) );
        ann.setSource( new String( getSource() ) );
        if ( getProperties() != null ) {
            ann.setProperties( ( PropertiesMap ) getProperties().copy() );
        }
        else {
            ann.setProperties( null );
        }
        return ann;
    }

    public Confidence getConfidence() {
        return _confidence;
    }

    public String getDesc() {
        return _desc;
    }

    public String getEvidence() {
        return _evidence;
    }

    public PropertiesMap getProperties() {
        return _properties;
    }

    public String getRef() {
        return _ref;
    }

    public String getSource() {
        return _source;
    }

    public String getType() {
        return _type;
    }

    public Uri getUri() {
        return _uri;
    }

    private void init() {
        _desc = "";
        _type = "";
        _source = "";
        _ref = "";
        _evidence = "";
        _confidence = null;
        _properties = null;
        _uri = null;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( ForesterUtil.isEmpty( getDesc() ) ) {
            return false;
        }
        return getDesc().equals( ( ( Annotation ) data ).getDesc() );
    }

    public void setConfidence( final Confidence confidence ) {
        _confidence = confidence;
    }

    public void setDesc( final String desc ) {
        _desc = desc;
    }

    public void setEvidence( final String evidence ) {
        _evidence = evidence;
    }

    public void setProperties( final PropertiesMap property ) {
        _properties = property;
    }

    public void setRef( final String ref ) {
        _ref = ref;
    }

    public void setSource( final String source ) {
        _source = source;
    }

    public void setType( final String type ) {
        _type = type;
    }

    public void setUri( final Uri uri ) {
        _uri = uri;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( NHXtags.ANNOTATION ); //TODO change on NHXv2!!!!! //TODO
        sb.append( ForesterUtil.replaceIllegalNhxCharacters( getDesc().toString() ) );
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( ( getConfidence() != null ) || ( getProperties() != null ) || ( getUri() != null )
                || !ForesterUtil.isEmpty( getDesc() ) ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( indentation );
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.ANNOTATION,
                                          PhyloXmlMapping.ANNOTATION_REF_ATTR,
                                          getRef(),
                                          PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR,
                                          getEvidence(),
                                          PhyloXmlMapping.ANNOTATION_TYPE_ATTR,
                                          getType(),
                                          PhyloXmlMapping.ANNOTATION_SOURCE_ATTR,
                                          getSource() );
            if ( !ForesterUtil.isEmpty( getDesc() ) ) {
                PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.ANNOTATION_DESC, getDesc(), indentation );
            }
            if ( getConfidence() != null ) {
                getConfidence().toPhyloXML( writer, level, indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
            }
            if ( getProperties() != null ) {
                getProperties().toPhyloXML( writer, level, indentation );
            }
            if ( getUri() != null ) {
                getUri().toPhyloXML( writer, level, indentation );
            }
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( indentation );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.ANNOTATION );
        }
        else {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.ANNOTATION,
                                             PhyloXmlMapping.ANNOTATION_REF_ATTR,
                                             getRef(),
                                             PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR,
                                             getEvidence(),
                                             PhyloXmlMapping.ANNOTATION_TYPE_ATTR,
                                             getType(),
                                             PhyloXmlMapping.ANNOTATION_SOURCE_ATTR,
                                             getSource(),
                                             indentation );
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
