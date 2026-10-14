# Channel Estimator

Repositorio correspondiente al desarrollo del **Estimador de Canal Adaptativo** dentro del proyecto **FI-AFE**.

Este repositorio se utilizará para el desarrollo, simulación y validación del estimador de canal.

El trabajo se organizará mediante **GitHub Projects, Issues, branches y Pull Requests**, de forma que cada tarea tenga un responsable y se desarrolle de manera independiente antes de integrarse al proyecto.

## Documentación

La documentación general del proyecto FI-AFE se encuentra disponible en la Wiki del repositorio principal:

[FI-AFE_2026 Wiki](https://github.com/Fundacion-Fulgor/FI-AFE_2026/wiki)

La documentación específica del estimador de canal se mantendrá en este repositorio.

## Estructura inicial

El repositorio se organiza por módulos, agrupando en cada uno sus modelos, simulaciones y archivos de implementación.

En esta primera etapa, el desarrollo estará enfocado en la simulación de alto nivel en Python.

```text
channel_estimator/
│
├── top/
│   ├── python/
│   └── simulations/
│
├── modules/
│   ├── signal_generation/
│   │   ├── python/
│   │   └── simulations/
│   │
│   ├── analog_channel/
│   │   ├── python/
│   │   └── simulations/
│   │
│   ├── cal_adc/
│   │   ├── python/
│   │   └── simulations/
│   │
│   ├── channel_estimator/
│   │   ├── python/
│   │   └── simulations/
│   │
│   └── correlator/
│       ├── python/
│       └── simulations/
│
├── docs/
├── .gitignore
├── requirements.txt
└── README.md
```

### Descripción de las carpetas

- `top/`: integración de los distintos módulos y simulaciones del sistema completo.
- `modules/`: desarrollo independiente de los módulos que componen el sistema y su entorno de validación.
- `docs/`: documentación técnica, decisiones de diseño y resultados relevantes.

Dentro de `modules/` se encuentran:

- `signal_generation/`: generación de señales de prueba, incluyendo PRBS, mapper PAM4 y pulse shaping.
- `analog_channel/`: modelo del canal analógico equivalente SISO PTV, con 16 kernels FIR.
- `cal_adc/`: modelo comportamental del ADC de calibración y sus no idealidades.
- `channel_estimator/`: desarrollo del estimador adaptativo, incluyendo selección de fases, ventanas de referencia y LMS.
- `correlator/`: modelo para el alineamiento temporal entre la referencia digital y las observaciones del Cal ADC.

Cada módulo contiene:

- `python/`: implementación de los modelos de alto nivel.
- `simulations/`: simulaciones, pruebas funcionales y validación del módulo.

La carpeta `top/` permite integrar los modelos y realizar simulaciones del sistema completo.

Los módulos de generación de señales, canal analógico y Cal ADC forman parte del entorno de simulación necesario para desarrollar y validar el estimador.

La estructura se ampliará a medida que avance el proyecto, incorporando modelos en punto fijo, RTL y verificación dentro de los módulos correspondientes.
## Organización del trabajo

Las tareas se organizarán utilizando **GitHub Projects**.

Para cada tarea:

1. Se creará un Issue y se asignará al colaborador correspondiente.
2. El colaborador creará una branch asociada al Issue.
3. La nueva branch deberá partir de `develop`.
4. Todo el desarrollo de la tarea deberá realizarse en esa branch.
5. Los commits deberán ser descriptivos.
6. Al finalizar, se deberá crear un Pull Request hacia `develop`.
7. El Pull Request será revisado antes de integrarse.
8. Si se solicitan cambios, deberán realizarse sobre la misma branch.
9. El Pull Request solo se integrará cuando esté aprobado.
10. Luego del merge, se cerrará el Issue y se eliminará la branch.

El flujo general será:

```text
GitHub Project
      ↓
Issue asignado
      ↓
Branch desde develop
      ↓
Desarrollo y pruebas
      ↓
Pull Request hacia develop
      ↓
Review
      ↓
Merge
      ↓
Cierre del Issue
      ↓
Eliminación de la branch
```

Las branches representan tareas específicas.

Ejemplo:

```text
Issue #12
Implement calibration ADC model

            ↓

12-implement-calibration-adc-model

            ↓

Pull Request → develop

            ↓

Review y aprobación
```

## Seguimiento de tareas

Cada Issue deberá incluir:

- Objetivo de la tarea.
- Actividades previstas.
- Responsable asignado.
- Avances realizados.
- Pruebas ejecutadas y resultados obtenidos.
- Dificultades o decisiones pendientes.

A medida que avance el desarrollo, se deberán actualizar las actividades completadas y registrar los resultados relevantes.

Los resultados de las simulaciones, gráficos y pruebas deberán estar disponibles en el repositorio o vinculados desde el Issue o Pull Request correspondiente.

El objetivo es mantener un registro de **qué se implementó, qué se validó y qué queda pendiente**.

## Branches principales

El repositorio utiliza dos branches principales:

### `main`

Contiene las versiones estables y validadas del proyecto.

No se debe desarrollar ni realizar push directamente sobre esta branch.

### `develop`

Es la branch de integración del proyecto.

Todas las branches de trabajo deben crearse a partir de `develop` y los cambios deben volver a `develop` mediante Pull Request.

No se debe desarrollar ni realizar push directamente sobre esta branch.

## Pull Requests

Todos los cambios que se quieran integrar a `develop` deberán realizarse mediante Pull Request.

Los colaboradores no deben mergear sus propios Pull Requests.

Cada Pull Request debe indicar brevemente:

- Qué tarea resuelve.
- Qué cambios se realizaron.
- Cómo fueron probados.
- Qué Issue está asociado.
- Si quedó algo pendiente.

Para relacionar el Pull Request con el Issue correspondiente, incluir en su descripción:

```text
Refs #12
```

Reemplazando `12` por el número del Issue.

Una vez aprobado e integrado el Pull Request en `develop`, se verificará que la tarea esté completa y se cerrará el Issue.

**Nota:** GitHub realiza el cierre automático mediante palabras clave como `Closes #12` cuando los cambios se integran en la rama predeterminada del repositorio. Como el flujo de trabajo utiliza `develop` como rama de integración, los Issues se cerrarán manualmente después del merge a `develop`, salvo que se configure otro mecanismo.

## Proyecto FI-AFE

Este repositorio forma parte del proyecto general:

[FI-AFE_2026](https://github.com/Fundacion-Fulgor/FI-AFE_2026)

El módulo se integra al repositorio principal como un submódulo del sistema completo.
