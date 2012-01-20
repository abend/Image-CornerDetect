

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'Image::CornerDetect',
    VERSION_FROM => 'lib/Image/CornerDetect.pm', # finds \$VERSION
    AUTHOR       => 'Sasha Kovar (sasha-cpan@arcocene.org)',
    ABSTRACT     => 'Detect corners in images via the Xiaochen He method.',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
