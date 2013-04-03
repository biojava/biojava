// $Id: SurfacingConstants.java,v 1.5 2008/12/20 01:42:02 cmzmasek Exp $
//
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

package org.forester.surfacing;

import org.forester.util.ForesterUtil;

public class SurfacingConstants {

    public static final String GOOGLE_WEB_SEARCH_LINK       = "http://www.google.com/search?q=";
    public static final String GOOGLE_SCHOLAR_LINK          = "http://scholar.google.com/scholar?q=";
    public static final String GOOGLE_SCHOLAR_LIMITS        = "&as_subj=bio&as_subj=med&as_subj=chm&num=100";
    public static final String AMIGO_LINK                   = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    public static final String PFAM_FAMILY_ID_LINK          = "http://pfam.sanger.ac.uk/family?id=";
    public static final String NL                           = ForesterUtil.LINE_SEPARATOR;
    public static final String TAXONOMY_LINK                = "http://beta.uniprot.org/taxonomy/?query=";
    static final boolean       SECONDARY_FEATURES_ARE_SCOP  = true;
    static final String        SECONDARY_FEATURES_SCOP_LINK = "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?key=";
    public static final String NONE                         = "[none]";
    public static final String UNIPROT_LINK                 = "http://beta.uniprot.org/taxonomy/?query=";
    public static final String GO_LINK                      = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    public static final String EOL_LINK                     = "http://www.eol.org/search?q=";
    public static final String TOL_LINK                     = "http://www.googlesyndicatedsearch.com/u/TreeofLife?q=";
    public static final String WIKIPEDIA_LINK               = "http://wikipedia.org/wiki/";
}
