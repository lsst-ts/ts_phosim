#!/usr/bin/env groovy

pipeline {

    agent {
        // Use the docker to assign the Python version.
        // Use the label to assign the node to run the test.
        // It is recommended by SQUARE team do not add the label and let the
        // system decide.
        docker {
            image 'lsstts/develop-env:develop'
            args '--entrypoint="" -u root'
        }
    }

    options {
      disableConcurrentBuilds()
      skipDefaultCheckout()
    }

    triggers {
        pollSCM('H * * * *')
    }

    environment {
        // Position of LSST stack directory
        LSST_STACK = "/opt/lsst/software/stack"
        // Pipeline stack version
        STACK_VERSION = "current"
        // XML report path
        XML_REPORT = "jenkinsReport/report.xml"
        // Module name used in the pytest coverage analysis
        MODULE_NAME = "lsst.ts.phosim"
        // PlantUML url
        PLANTUML_URL = "https://managedway.dl.sourceforge.net/project/plantuml/plantuml.jar"
        // Authority to publish the document online
        user_ci = credentials('lsst-io')
        LTD_USERNAME = "${user_ci_USR}"
        LTD_PASSWORD = "${user_ci_PSW}"
        DOCUMENT_NAME = "ts-phosim"
    }

    stages {

        stage ('Cloning Repos') {
            steps {
                dir(env.WORKSPACE + '/ts_phosim') {
                    checkout scm
                }
                dir(env.WORKSPACE + '/phosim_utils') {
                    git branch: 'master', url: 'https://github.com/lsst-dm/phosim_utils.git'
                }
                dir(env.WORKSPACE + '/ts_wep') {
                    git branch: 'master', url: 'https://github.com/lsst-ts/ts_wep.git'
                }
                dir(env.WORKSPACE + '/ts_ofc') {
                    git branch: 'master', url: 'https://github.com/lsst-ts/ts_ofc.git'
                }
            }
        }

        stage ('Building the Dependencies') {
            steps {
                // When using the docker container, we need to change
                // the HOME path to WORKSPACE to have the authority
                // to install the packages.
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source ${env.LSST_STACK}/loadLSST.bash

                        cd phosim_utils/
                        setup -k -r . -t ${env.STACK_VERSION}
                        scons

                        cd ../ts_wep/
                        setup -k -r .
                        scons python
                    """
                }
            }
        }

        stage ('Unit Tests and Coverage Analysis') {
            steps {
                // Direct the HOME to WORKSPACE for pip to get the
                // installed library.
                // 'PATH' can only be updated in a single shell block.
                // We can not update PATH in 'environment' block.
                // Pytest needs to export the junit report.
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source ${env.LSST_STACK}/loadLSST.bash

                        cd phosim_utils/
                        setup -k -r . -t ${env.STACK_VERSION}

                        cd ../ts_wep/
                        setup -k -r .

                        cd ../ts_ofc/
                        setup -k -r .

                        cd ../ts_phosim/
                        setup -k -r .
                        pytest --cov-report html --cov=${env.MODULE_NAME} --junitxml=${env.WORKSPACE}/${env.XML_REPORT}
                    """
                }
            }
        }
    }

    post {
        always {
            // The path of xml needed by JUnit is relative to
            // the workspace.
            junit "${env.XML_REPORT}"

            // Publish the HTML report
            publishHTML (target: [
                allowMissing: false,
                alwaysLinkToLastBuild: false,
                keepAll: true,
                reportDir: 'ts_phosim/htmlcov',
                reportFiles: 'index.html',
                reportName: "Coverage Report"
            ])

            script{
              withEnv(["HOME=${env.WORKSPACE}"]) {
                def RESULT = sh returnStatus: true, script: """
                  source ${env.LSST_STACK}/loadLSST.bash

                  curl -O ${env.PLANTUML_URL}

                  pip install sphinxcontrib-plantuml

                  source ${env.LSST_STACK}/loadLSST.bash

                  cd phosim_utils/
                  setup -k -r . -t ${env.STACK_VERSION}

                  cd ../ts_wep/
                  setup -k -r .

                  cd ../ts_ofc/
                  setup -k -r .

                  cd ../ts_phosim/
                  setup -k -r .

                  package-docs build

                  pip install ltd-conveyor

                  ltd upload --product ${env.DOCUMENT_NAME} --git-ref ${env.GIT_BRANCH} --dir doc/_build/html
                    """

              if ( RESULT != 0 ) {
                  unstable("Failed to push documentation.")
              }
            }
          }
        }

        cleanup {
            // Change the ownership of workspace to Jenkins for the clean up
            // This is to work around the condition that the user ID of jenkins
            // is 1003 on TSSW Jenkins instance. In this post stage, it is the
            // jenkins to do the following clean up instead of the root in the
            // docker container.
            withEnv(["HOME=${env.WORKSPACE}"]) {
                sh """
                    chown -R 1003:1003 ${HOME}/
                """
            }
            // clean up the workspace
            deleteDir()
        }
    }
}
